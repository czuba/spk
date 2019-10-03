function [spk, OpenedFileName] = plx_readerPar_Pldaps(OpenedFileName, arg, sync)
% Reads in the data from PLX & PL2 files using parallel toolbox. 
% based on the SDK provided by Plexon
% 
% --------------------------
% 2018-05 update notes, TBC:
%   Lots of this may be obsolete now that most all sorting is done by
%   kilosort, which have a completely distinct output structure/interaction.
%   But still lots of good code practices, useful data structures, and
%   parallel processing examples live in here that should be carried over
%   to PLDAPS data infrastructure.
% --------------------------
% 
%   [OpenedFileName]    Full path to .plx/.pl2 file (...on Mac, will expand "~" to home dir)
%   [arg] flag options  (...hackily complex, but is what it is)
%       default: [666,0]  Loads all spikes & collapse sort codes, no waveforms;
% 
%       arg(1) == 0:    Loads all spikes (including unsorted), retain existing sortcodes.
%       arg(1) == 1:    Ignore unsorted spikes (PLX sortcode == 0)
%       arg(1) == 666:  Collapse all spike sort codes
%                           (good for ignoring junk units created during recording)
%       arg(2) == 0:    Dont get waveforms.
%       arg(2) == 1:    Load waveforms, but only pass out wf mean & ci of each unit
%       arg(2) == 2:    Load waveforms, pass out wf mean, ci, & every individual wf (this will be big!)
%       arg(2) == 3;    Load waveforms, pass out wf mean & means from 4 temporal intervals across file
%                           (good for coarse stability check)
% 
% 
% [spk] data structure formatted as:  (see plx_readerPar)
%           id:         unit id (chan*10 + sortUnit#)
%           ts:         timestamp cell (nUnits x 1; in ms [if sr==1000])
%           sync:       sync times  (in ms [if sr==1000])
%           snr:        waveform snr (nUnits x 1)
%           wf:         struct with mean and quantiles of waveforms (fields mu & ci)
%                       * can also request to have .all field with all waveforms...will make file huge & not recommended
%           trodality:  trodal recording mode reported by PLX
%           info:       any comment text entered into Plexon at time of recording
%           n:          nSpikes per unit (nUnits x 1)
% 
% Plot waveforms as: (where n=desired unit index)
%       plot(wf.mu(:,n),'-'); hold on, plot(wf.ci(:,:,n), '--')
% 
% NOTE: Currently only gets one event channel (Expo-centric). Pldaps provides
%       much more event info & chans that this needs to be updated for
% 
% 2013-xx-xx  Written by T.Czuba, based on Plexon SDK
% 2013-2017   TBC  Various edits, extensions, output struct developments
% 2017-10-13  TBC  Commented & cleaned up a bit.

if nargin <2 || isempty(arg)
    arg = [666, 0];
elseif length(arg)<2
    arg = [arg,0];
end

if nargin <1 || isempty(OpenedFileName)
    [fileName, filePath] = uigetfile('*.*', 'Select spike source file');
    if ~fileName
        % User pressed cancel
        spk = [];
        return
    end
    OpenedFileName = fullfile(filePath, fileName);
end

% Initialize sync with PLDAPS file
if nargin>=3 && isstruct(sync)
    dosync = true;
    plxInfo = sync.info;
else
    dosync = false;
    fprintf('\tGetting plx file info...\n')
    plxInfo = getPlxInfo(OpenedFileName);
    sync = [];
end

doalign = 1;
plxInfo.alignWfs = doalign;

% expand home directory path (else will crash in mexplex)
if strcmp(OpenedFileName(1),'~')
    [~,jnk] = system('echo $HOME');  jnk = jnk(1:end-1);
    OpenedFileName = fullfile(jnk,OpenedFileName(2:end));
end
    
% % % [OpenedFileName, Version, Freq, Comment, Trodalness, NPW, PreThresh, SpikePeakV, SpikeADResBits, SlowPeakV, SlowADResBits, Duration, DateTime] = plx_information(OpenedFileName);
% % % [tscounts, wfcounts, evcounts, slowcounts] = plx_info(OpenedFileName,1);
[OpenedFileName, ~, datFreq, comment, trodality, npw, ~, ~, ~, ~, ~, ~, ~] = plx_information(OpenedFileName);
[tscounts, ~, evcounts, ~] = plx_info(OpenedFileName,1);


% skip unsorted spikes if flagged
if arg(1)==1
    tscounts(1,:) = 0;
end

% get event chanmap
[~,evchans] = plx_event_chanmap(OpenedFileName);

% Determine trodality and spike indices accordingly
ss = find(sum(tscounts));
if trodality >1
    sel = ss(ismember(ss, 2:trodality:(length(tscounts))));
else
    sel = ss;
end
% trim redundancy (present if trodality>1)
tscountsn = zeros(size(tscounts));
tscountsn(:,sel) = tscounts(:,sel);
clear sel ss;

[xx,yy] = find(tscountsn);
nunits = length(xx);

% Initialize variables for output
[id, snr] = deal(nan(nunits,1));
[ts, wfmu, wfci] = deal(cell(nunits,1));

%waveform ci levels (outside of parfor loop!)
ciLvl = [.16 .84];
if arg(2)==3 % coopt nci for n temporal averaging windows
    nci = 8;
else
    nci = length(ciLvl);
end


%% Get the spikes!
fprintf('\tLoading spike times (and wfs)...\n')

% start parallel pool if not already present
if isempty(gcp('nocreate'))
    parpool;
end

% progress in command window
fprintf(['\n[',repmat('.', [1,nunits]),']\n[ '])

tic
parfor i = 1:nunits
    id(i) = (yy(i)-1)*10 + (xx(i)-1);       % spike id = chan*10 + unit#
    id(i) = ceil((yy(i)-1)/trodality)*10 + (xx(i)-1);       % spike id = chan*10 + unit#

    if ~arg(2)                %#ok<PFBNS>   % just get timestamps
        [~, ts{i}] = plx_ts(OpenedFileName, yy(i)-1 , xx(i)-1 );
    else % get timestamps and waveforms
        if trodality==1
            [~, npw, ts{i}, wv] = plx_waves_v(OpenedFileName, yy(i)-1 , xx(i)-1 );
        else            
            wv = cell([1,trodality]);  % clear/setup wv variable to receive each chan of stereo/tetrode
            [~, npw, ts{i}, wv{1}] = plx_waves_v(OpenedFileName, yy(i)-1 , xx(i)-1 );
            for j = 2:trodality
                [~, ~, ~, wv{j}] = plx_waves_v(OpenedFileName, yy(i)-1 +(j-1), xx(i)-1 );
            end
            wv = cell2mat(wv);
        end

        % compute mean and quantiles of wf
        if size(wv,1)>1            
            % align waveforms (see sub-function)
            if doalign
                [wv] = alignWfs(wv, npw*trodality, 10);
            end
            
            wfmu{i} = nanmean(wv);
            if numel(ciLvl)<nci
                % Hack to extract evenly spaced quantiles
                % instead of typical +/-1 standard deviation (see "if arg(2)==3" below)
                wfci{i} = quantile(wv, nci);
            else
                wfci{i} = quantile(wv, ciLvl);    %[.025 .16 .84 .975]
            end
        else % don't crash if only one spike in unit
            wfmu{i} = single(wv);
            wfci{i} = single(nan(nci, size(wv,2)));
        end
        
        % make wfci the means of temporal quartiles if requested with arg(2)==3
        if arg(2)==3
            % get timepoints for windows (ok/intentional that not abutting)
            jnkt = ceil(linspace(1,length(ts{i}),nci+1));
            for ti = 1:nci
                wfci{i}(ti,:) = nanmean( wv(jnkt(ti):jnkt(ti+1), :));
            end
        end
        
        % only pass out whole pile of waveforms if requested with arg(2)==2
        if arg(2)==2
            wfall{i} = wv;
        end
        
        % compute snr
        A = max(wfmu{i}) - min(wfmu{i});
        wv = bsxfun(@minus, wv, wfmu{i});
        wv = std(wv(~isnan(wv)));
        snr(i) = A/(2*wv);
        
    end
    % print progress
    fprintf('\b|\n')
end
toc
fprintf('\b]\n%3.2f sec per unit\n\n',toc/nunits)

if dosync   %nargin>=3 && isstruct(sync)
    % Ensure any time conversions applied to syncs are also applied to spikes
    if (sync.plxTrialTs(1)-sync.info.pl2ptbFit(2))>0    %abs(sync.plxTrialTs(1)-dv.pds.data{1}.trstart) < abs(dv.uprb.sync.info.pl2ptbFit(2))
        %   if spike syncs appear to have been converted during syncing (default), then convert spike timestamps to PTB time too
        ts = cellfun(@(x) sync.plx2ptb(x), ts, 'uni',0);
    end

    % trim to trial range +- 1 min
    % Sets range of spikes to include in output
    tsLim = [min(sync.strobe(:,1)), max(sync.strobe(:,1))] + [-60 60];
    for i = 1:size(ts)
        ii = ts{i}<tsLim(1) | ts{i}>tsLim(2);
        ts{i}(ii) = [];
        if exist('wv','var')
            keyboard % ...this may still need updating for PLDAPS, but also may not be worth it given kilosort pipeline --TBC 2018-10
            wv{i}(ii) = [];
        end
    end
end

%% Syncs/events  ~~disabled. just deal with spikes here & passthrough sync~~
if 0
    % Needs to be updated to produce cell of sync event channels
    [~, sync] = plx_event_ts(OpenedFileName, evchans(evcounts==max(evcounts)));
    if sum(sync)<=0
        % evchan mislabeled...warn and seek alternates
        fprintf('\t!Error in evchannel mapping!\t')
        
        if strcmp(OpenedFileName(end-2:end), 'pl2')
            % dig into PL2 struct b/c OFS V3 likes to garble eventTS after saving/exporting
            fprintf('(Likely PL2 headers got screwed up by perpetual bugs in Offline Sorter V3)\n\t...attempting to resolve...')
            pl2 = PL2GetFileIndex(OpenedFileName);
            evcounts = cellfun(@(x) x.NumEvents, pl2.EventChannels)';
            
            ech = find(evcounts==max(evcounts));
            jnk = PL2EventTs(pl2.FilePath, ech(1));
            sync = jnk.Ts; clear jnk
            
        else
            %else lumber through plx evchans till find them.
            fprintf('(Likely an artifact of combining plx files with different sync chans)\n\t...attempting to resolve...')
            for i = (max([1,min(evchans)-5])):(max(evchans)),
                evcounts(i) =  plx_event_ts(OpenedFileName, i);
            end
            ech = find(evcounts==max(evcounts));
            [~, sync] = plx_event_ts(OpenedFileName, ech(1));
        end
        
        if sum(sync)<=0
            warning('Zero syncs found...this is probably bad.')
        else
            warning('\n\tFound alt evchan# == %d when zero syncs were found with default event chan mapping. Be warned...\n', ech)
        end
    end
end
    
if arg(1) == 666
    % kill all units. Pile into unsorted chans.
    [u, ~, ui] = unique(floor(id./10));
    ts0 = cell(size(u));
    for i = 1:length(u)
        ts0{i} = sort(cell2mat(ts(ui==i)));
    end
    % Need to redo snr too
    % -- Fix size here, approximate snr from weighted sum below
    snr0 = nan(size(ts0));
    
    % collapse wfs too (by weighted means)
    if arg(2)
        for i = 1:length(u)
            wwf = cellfun(@numel, ts(ui==i));
            wwf = wwf/norm(wwf,1);
            if sum(ui==i)>1
                wfmu0{i} = sum(bsxfun(@times, cat(1, wfmu{ui==i}), wwf));
                jnk = bsxfun(@times, cat(1, wfci{ui==i}), kron(wwf, ones(nci, 1)));    %ones(numel(wwf)*nci/2, 1)));
                for ii = 1:nci%size(jnk,1)/2
                    wfci0{i}(ii,:) = sum(jnk(ii:2:end, :));
                end
            else
                wfmu0{i} = wfmu{ui==i};
                wfci0{i} = wfci{ui==i};
            end
            % approx hack of snr
            %   (surely invalid, but this mode is non-quantitative smash & grab)
            snr0(i) = sum(bsxfun(@times, snr(ui==i), wwf));

        end
        wfmu = wfmu0'; % match cell dims to orig
        wfci = wfci0';
        snr = snr0;
    end
    
    id = u*10;  % keep id in same chan*10 format so can tell later that these are unsorted spikes
    ts = ts0;
    nunits = numel(ts);
    clear ts0 wfci0 wfmu0 snr0

end

% Package up id,ts,sync
spk.id    = id;
spk.ts    = ts;
% spike count for each unit
spk.n = cellfun(@length, spk.ts);

% Sync structure
spk.sync  = sync;

% if got snr/waveforms too
if arg(2)
    spk.snr   = snr;
    spk.wf.mu = cell2mat(wfmu)';
    % order wf confidence intervals (act. 1 & 2 sd quantiles) so they will plot like wf.mu
    %       i.e.   plot(wf.mu(:,n),'-'); hold on, plot(wf.ci(:,:,n), '--')
    spk.wf.ci = permute( reshape(cell2mat(wfci), [nci, nunits, npw*trodality]), [3,1,2]);
        
    % may request all wfs too, but not a default output (too big)
    if arg(2)==2
        spk.wf.all = wfall;
    end
    
else % empty placeholders
    spk.snr   = snr;
    spk.wf.mu = cell2mat(wfmu)';
    spk.wf.ci = cell2mat(wfci)';     
end

spk.trodality = trodality;
spk.info = plxInfo;

end

% % % % % % % % % % % % 
% % % SubFunctions
% % % % % % % % % % % % 

%% alignWfs
function [wf2, lag] = alignWfs(wf, wEnd, maxShift)
% function [wf2, lag] = alignWfs(wf, wEnd, maxShift)
% 
% Aligns waveforms by estimating the best fitting delay/shift (crosscorrelation via finddelay.m)
%   [wf]        = Matrix of waveforms, as returned by plx_waves
%   [wEnd]      = Meant to clip wf length
%                   ~!~ NOT implimented ~!~  Defaults to all points in wf: size(wf,2)   TBC 2014-12-01
%                   --minimize noise/var from wf length exceeding spike width (e.g. 1200µs) and clipped wf tails
%                   --ideally would really be range of wf points to use, but implementation is shifty b/c of device
%                     dependent sampling rate...
%   [maxShift]  = Maximum waveform shift (in wf samples [unfortunately] NOT ms)
%                   --Default == 10
% 
% 2014-09-XX TBC  Wrote it
% 2014-12-01 TBC  Commented

[nwf, nsamp] = size(wf);

if ~exist('wEnd','var') || isempty(wEnd)
    wEnd = nsamp;
end
if ~exist('maxShift','var') || isempty(maxShift)
    maxShift = 10;
end

% find align shifts
wfmu = nanmean(wf);
lag = -finddelay(wf', repmat(wfmu,[nwf,1])', maxShift)';

% init nan array (singles save memory)
wf2 = nan(size(wf), 'single');

% do the shifting
for i = 1:nwf,
    if lag(i)>=0,
        wf2(i, 1:nsamp-lag(i)) = wf(i, lag(i)+1:nsamp);
    else
        wf2(i, 1-lag(i):nsamp) = wf(i, 1:(nsamp+lag(i)));
    end

end
end
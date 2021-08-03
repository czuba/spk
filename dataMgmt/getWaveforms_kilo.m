function [wf, dupes, wfFigs] = getWaveforms_kilo(tsIdx, varargin)
% function wf = getWaveforms_kilo(tsIdx, pars)
%
% Extracts individual spike waveforms from the raw datafile, for multiple
% clusters. Returns the waveforms and their means within clusters.
%
% based on getWaveforms from https://github.com/cortex-lab/spikes.
% ...extended & merged with Czuba's Einstein pipeline
%
% % EXAMPLE INPUT
% tsIdx = cell of spike timestamps in samples;  {nUnits, 1}(nspikes,1))
% 
% pars.dataDir = '/path/to/data/';    % KiloSort/Phy output folder
% pars.fileName = 'data.dat';         % .dat file containing the raw voltage data (this should be BP filtered)
% pars.dataType = 'int16';            % Data type of .dat file (always int16)
% pars.nCh = 32;                      % Number of channels that were streamed to disk in .dat file
% pars.wfWin = [-40 41];              % Number of samples before and after spiketime to include in waveform
% pars.nWf = 2000;                    % Number of waveforms per unit to use for averaging
% pars.align = true;                  % Improve waveform alignment using /spk/dataMgmt/alignWfs_kilo.m
%                                     % (coarse trough alignment, then fine crosscorrel. align to mean)
%
% % OUTPUT
% wf.mu         % [nSamples, nCh, nUnits]       Average of all waveforms (per channel)
% wf.ci         % [nSamples, nCh, nCi, nUnits]  Conf intervals ( [.16,.84] prctiles)
% wf.peakCh     % channel # of peak waveform amplitude for each unit
% wf.snr        % Waveform SNR (only computed if pars.align==true[default]; invokes  /spk/dataMgmt/getWaveforms_kilo.m)
% 
% wf.pars       % parameters used to constrcut [wf] output
%
% nUnits:   number of different units/clusters  == size(tsIdx,1) == numel(tsIdx) 
% nSamples: number of samples per waveform      == diff(pars.wfWin)
% 
% % USAGE
% wf = getWaveForms_kilo(tsIdx, varargin);


%% Parse inputs & setup default parameters
pp = inputParser();
% Basic data parameters
pp.addParameter('dataDir', pwd, @isfolder)
pp.addParameter('fileName', []);
pp.addParameter('dataType', 'int16');
pp.addParameter('raw2v', []);
pp.addParameter('nCh', []);         % Number of channels that were streamed to disk in .dat file
pp.addParameter('trodality', 1);    % if >1, will reshape waveforms to n*wfSamples
pp.addParameter('chDepth', []);     % depth of channel; input from spk.info.chanMap.ycoords (def 1:nCh)
pp.addParameter('sr', 40000);       % Raw data sampling rate (def 40,000)
pp.addParameter('nWf', 2000);       % nWaveforms to load/average
pp.addParameter('wfWin', [-500, 1000]*1e-6);  % waveform samples [before, after] timestamp (in seconds OR samples; should be 1:2 ratio relative to 0)
pp.addParameter('ci', [.16 .84]);   % waveform ci levels (def +/-1sd)
pp.addParameter('all', false);      % return all waveforms loaded (i.e. not just mean & ci wfs)     % Soon to be depricated
pp.addParameter('align', true);     % align waveforms (used for SNR calc & waveform mean/ci plots)
pp.addParameter('acg', true);       % compute autocorrelogram (uses CCG.m/mex for speed, from FMAToolbox)
pp.addParameter('plotSnr', true);   % 


pp.parse(varargin{:});

pars = pp.Results;

% Fill in missing inputs
if isempty(pars.fileName)
    % load info from params.py
    kiloInfo = loadParamsPy( fullfile(pars.dataDir, 'params.py'));
    pars.fileName = kiloInfo.dat_path;  % .dat file containing raw voltages (can be relative to .dataDir; e.g. '../../raw/file.dat')
    pars.dataType = kiloInfo.dtype;     % Data type of .dat file (this should be BP filtered)
    pars.nCh = kiloInfo.n_channels_dat; % Number of channels that were streamed to disk in .dat file
    pars.sr = kiloInfo.sample_rate;
end

pars.datPath = fullfile(pars.dataDir, pars.fileName);
if ~exist(pars.datPath,'file')
    % try to prepend raw file location in standard hierarchy (workaround for symlinks being disallowed on certain file systems)
    pars.datPath = fullfile(pars.dataDir, '..','..','raw', pars.fileName);
end

% Convert wf duration to samples
if abs(pars.wfWin(1))<1 || all(round(pars.wfWin)~=pars.wfWin)
    pars.wfWin = round(pars.wfWin*pars.sr);
end
if isempty(pars.chDepth), pars.chDepth = 1:pars.nCh; end

% create output struct
wf = struct('mu',[], 'ci',[], 'pars', pars);

if ~exist(pars.datPath, 'file')
    % Exit gracefully if no raw dat file found
    fprintf(2, '~~~\tCould not find raw dat path to load waveforms from...continuing without them\n~~~\t%s\n', pars.datPath)
    return
end


%% Load .dat and KiloSort/Phy output
filenamestruct = dir(pars.datPath);
dataTypeNBytes = numel(typecast(cast(0, pars.dataType), 'uint8')); % determine number of bytes per sample
nSampTotal = filenamestruct.bytes/(pars.nCh*dataTypeNBytes);  % Number of samples per channel
datSize = [pars.nCh nSampTotal];
nSamples = length(pars.wfWin(1):pars.wfWin(end));

% Memory map raw dat file
mmf = memmapfile(pars.datPath, 'Format', {pars.dataType, datSize, 'x'});
% Confirm order in which data was streamed to disk; must be 1-indexed for Matlab
chMap = readNPY(fullfile(pars.dataDir, 'channel_map.npy'))+1;
nChInMap = numel(chMap);
nUnits = numel(tsIdx);


%% ACG parameters
binsz = .0002;
binDur = .015;
sr = pars.sr; % data sample rate (default==40k Hz)


%% Prepare vars for loop
% Initialize outputs
[allWfs, tmpMu, tmpCi, tmpPeakCh, wfSnr, wfFigs, acg, dupes] = deal(cell(nUnits,1));
[wfSnr{:}] = deal(nan); % initialize snr as nan
% Slice variables for parfor
nWf = pars.nWf;
nCh = pars.nCh;
raw2v = pars.raw2v;
wfWin = [pars.wfWin(1), pars.wfWin(end)];
ciLvls = pars.ci;
% Processing Flags
getAllWfs = pars.all; % Only compile ALL waveforms if requested (big & slow!)
doAlign = pars.align;
doAcg = pars.acg;
plotSnr = pars.plotSnr;

tt0 = tic;
fprintf('\nTotal units:\t[%s]\nLoading wfs:\t[', repmat('.',[nUnits,1]))

%% [par]for loop through each unit
% NOTE:  Can't use parfor if plotting
for u = 1:nUnits
    fprintf('.')
    tsU = tsIdx{u};
    theseWf = min([nWf, length(tsU)]);
    % select random subset of all timestamps
    wfi = randperm(length(tsU), theseWf);
    wf0 = tsU(wfi);
    
    % get waveform snippets from all channels
    tmpWf = getMemMapWfs(mmf, wf0, wfWin, 1:nCh);

    % find peak channel
    wfAmp = squeeze(range(nanmean(tmpWf,3)));
    [wfAmpMax, tmpPeakCh{u}] = max(wfAmp);
    
    if doAlign
        % compute alignment on ALL spikes in unit (crucial for detecting double-counted spikes)
        % - pass wfs to alignment function as:  [count, wfSamples]
        [~, lag] = alignWfs_kilo(squeeze(getMemMapWfs(mmf, tsU, wfWin, tmpPeakCh{u}))');
        % apply offset to all timestamps
        tsU = tsU+lag;  % % % PAUSE HERE TO CHECK FOR DOUBLE COUNTING AFTER ALIGNMENT % % %
        % resample [wf0] timestamps from aligned originals
        wf0 = tsU(wfi);
        % record logical index of double counts (defined as isi violations)
        dupes{u} = abs(diff(tsU))<=(0.0005*sr);
        % remove doubles (else will show up in acg)
        tsU(dupes{u}) = [];
    end
    % Retrieve [aligned] waveforms from ALL channels
    tmpWf = getMemMapWfs(mmf, wf0, wfWin, 1:nCh);
    % confirm peak channel (shouldn't change, but ensure consistent with output)
    wfAmp = squeeze(range(nanmean(tmpWf,3)));
    [wfAmpMax, tmpPeakCh{u}] = max(wfAmp);

    % Compute mean & conf. intervals on waveform
    tmpMu{u} = nanmean(tmpWf, 3);
    tmpCi{u} = quantile(tmpWf, ciLvls, 3);
    
    if doAlign
        %% compute waveform SNR
        wfSnr{u} = calcWaveformSNR(mmf, wf0, wfWin, 1:nCh, plotSnr); %tmpPeakCh{u}
        wfFigs{u} = gcf;
        set(wfFigs{u}, 'name',sprintf('u%03d snr%.2f',u,wfSnr{u}));
    end
    
    if doAcg
        %         % calculate acg while we're here (WITH alignment might be slowwww)
        %         [~, lag] = alignWfs_kilo(squeeze(getMemMapWfs(mmf, tsU, wfWin, tmpPeakCh{u}))');
        %         [acg{u}, acgX] = CCG((tsU+lag)/sr, ones(size(tsU)), 'binSize',binsz, 'duration',binDur);
        [acg{u}, acgX] = CCG((tsU)/sr, ones(size(tsU)), 'binSize',binsz, 'duration',binDur);
    end

    if getAllWfs
        allWfs{u} = tmpWf;
    end
end
fprintf(']\t Compiling outputs...')

% Expand parfor cell outputs
sz = size(tmpMu{1});
sz = [sz(1)*pars.trodality, sz(2)/pars.trodality];
wf.mu          = nan( [sz, nUnits]);
wf.ci          = nan( [sz, numel(pars.ci), nUnits]);
for u = 1:nUnits
    wf.mu(:,:,u) = reshape(tmpMu{u}, sz);
    wf.ci(:,:,:,u) = reshape(tmpCi{u}, [sz, numel(pars.ci)]);
end
wf.peakCh = cell2mat(tmpPeakCh);

if doAlign
    wf.snr = cell2mat(wfSnr);
    % record isi violation rate
    wf.isiErrRate = cellfun(@(x,d) sum(d)/numel(x), tsIdx, dupes);
end

if doAcg
    wf.acg.counts = cell2mat(acg');
    wf.acg.bins = acgX;
end

% Only deal with all wfs if we have to...
% Soon to be depricated (...far too big with typ. channel & spike counts, circa 2020)
if getAllWfs
    wf.all    = nan(nSamples, pars.nCh, pars.nWf, nUnits);
    for u = 1:nUnits
        wf.all(:,:,1:size(allWfs{u},3),u) = allWfs{u};
    end
end


% add mean waveform [& ACG] plots
if plotSnr && doAlign
    for u = 1:nUnits
        % select figure
        figure(wfFigs{u})
        spx = 2;
        spy = 3;
        % ax = findobj(wfFigs{u},'type','axes');
        % figure(ax(1).Parent);
        % addWfInset(ax(1), wf, u, [0,-.5]);
        if ~doAcg
            addWfInset([spy,spx,3], wf, u, [], raw2v)
        else
            % plot mean waveform w/ci
            addWfInset([spy,spx,3], wf, u, [], raw2v)
            % plot ACG
            sp = subplot(spy,spx,5);
            plot(sp, wf.acg.bins, wf.acg.counts(:,u), '-','color',.5*[1 1 1])
            fmtPlot(sp,'small'); %box off, axis tight, yl = ylim; ylim(yl+.05*diff(yl)*[-1,1]), set(sp, 'fontsize',8);
            hold on;
            % reference lines
            plot(sp, [0,0;wf.acg.bins([1,end])']', [-100,10000;0,0]',':','color',.68*[1 1 1], 'xliminclude','off', 'yliminclude','off');
            ylabel('count')
            xlabel('autocorrelogram')
        end
        % Firing rate histogram
        sp = subplot(spy,spx,6);
        % this should also set bin limits == range of experiment, in case no spiking on ends
        histogram(tsIdx{u}(~dupes{u})/sr, 'Normalization','countdensity', 'DisplayStyle','stairs');
        fmtPlot(sp,'small');
        ylim(get(sp,'ylim').*[1,1.3]); % center histogram 
        %title('Firing rate')
        ylabel('spk/s')
        xlabel('expt time (sec)')
    end
end

fprintf('\n\t\tdone. (%3.1f sec)\n', toc(tt0))

end %main function


%% WF plotting example:
% 
%     figure;
%     nUnits = size(wf.mu,3);
%     spx = ceil(sqrt(nUnits));
%     showCi = 0;
% 
%     for u = 1:nUnits
%         h = subplot(spx, spx, u);
%         plot(h, wf.mu(:,:,u)); 
%         if showCi
%             hold on
%             jnk = wf.ci(:,:,:,u);
%             plot(h, jnk(:,:), '--');
%         end
%         box off
%         title(sprintf('unit %03d',u), 'fontsize',10)
%     end

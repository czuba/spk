function spk = getSpikes_kilo(kilopath, sync)
% 'dv' is standard data struct
%   -- Must have:
%       .paths.kilo field containing (full) path to spike datafile
%       -- can be either .nev or .plx (if sorted)
%       .expo field with standard expo struct
%   -- Returns resulting spike and sync timestamps in spk struct
% 
%   flag ==
%       0: get spikes only
%       1: get spikes, calc snr     [default]
%       2: calc snr only
% 

% defaults
if ~exist('kilopath','var') || isempty(kilopath)
    kilopath = [];
end

% Initialize sync with PLDAPS file
if nargin>=1 && isstruct(sync)
    dosync = true;
    plxInfo = sync.info;
    % Sets range of spikes to include in output
    plxTrialTsSec = sync.ptb2plx(sync.plxTrialTs); %convert PLDAPS synced trial times to acquisition system time ('ptb'==psychtoolbox, 'plx'==plexon)
    tsSyncLimSecs = [mmin(plxTrialTsSec), mmax(plxTrialTsSec)] + [-1 120];
else
    dosync = false;
    %     fprintf('\tGetting plx file info...\n')
    plxInfo = []; % getPlxInfo(OpenedFileName);
    sync = [];
end


%% Kilosort info
% initial params
kiloInfo = loadParamsPy( fullfile(kilopath, 'params.py'));
% channel map (incl spatial coords of recording device)
kiloInfo.chanMap = load(fullfile(kilopath, 'chanMap.mat'));
% raw dat file info
kiloInfo.rawFile = fullfile(kilopath, kiloInfo.dat_path);   % TODO: use as destination for preprocessed [wf] struct
if ~exist(kiloInfo.rawFile, 'file')
    % try to prepend raw file location in standard hierarchy
    kiloInfo.rawFile = fullfile(kilopath, '..','..','raw',kiloInfo.dat_path);
end
% Load rawInfo file (created by czuba fork of Kilosort setup)
rawInfo = [kiloInfo.rawFile(1:end-4),'_rawInfo.mat'];
if exist(rawInfo, 'file')
    rawInfo = load(rawInfo);
else
    rawInfo = [];
end
kiloInfo.rawInfo = rawInfo;

% Carry over comment from plexon file (...redundancy in this, but not large enough to hurt)
if isfield(sync.info, 'comment')
    kiloInfo.comment = sync.info.comment;
else
    kiloInfo.comment = '';
end


%% Get Kilosort spike data
% load all spike times (in samples)
tsIdx = double(readNPY( fullfile(kilopath, 'spike_times.npy')));

% get template Id of each spike time
%   (this is a ZERO-based index to the cluster templates; **ADD 1** to make ONE-based!)
if exist(fullfile(kilopath, 'spike_clusters.npy'),'file')
    idSrc = 'spike_clusters.npy'; % only exists if manual curation performed
else
    idSrc = 'spike_templates.npy'; % initial unit outputs from kilosort
end
spikeId = double(readNPY( fullfile(kilopath, idSrc)));


%% Load Sort IDs (Nested Function)
% sortId == label given to cluster group during sorting
% Kilosort/Phy has default/special labels:  {'good','mua','unsorted','noise'},
% ...but users can potentially use any sorting name they choose.
% Split the difference here while keeping some consistency with Offline Sorter
%   0=unsorted, 1=good, 2=mua, 3+=userAlphabetical,  ...-1=noise?

spk.sortId = loadSortIds;% ({'noise','good'});


%% Format and trim

% Trim unwanted sort groups (def: 'noise')
%   -Do this here/early to potentially save a lot of time/space if excluded sortIds correspond to mmany 'spikes'
uid = [spk.sortId{:,3}];
good = ismember(spikeId, uid);
spikeId(~good) = [];
tsIdx(~good) = [];

% Adjust spikeId to a ONE-based index that we can use
[ksIds, ~, spikeId] = unique(spikeId);

% convert tsIdx & spikeIds into an nUnit-by-1 cell of timestamps (in seconds)
tsIdx = accumarray(spikeId, tsIdx, [], @(v) {v});    % tsIdx = {nUnits, 1}(nspikes,1))
% accumarray.m does not maintain order(see doc); sort outputs rather than risk misassignment of inputs
tsIdx = cellfun(@sort, tsIdx, 'uni',0);


%% Initial depth ordering (Hacky)
% Estimate depth of each unit based on mean waveform amplitude (adapted from github/cortex-lab/spikes, templatePositionsAmplitudes.m)
% ...inefficient, but better than waiting till later and needing to reorder mmany more fields derived from
% initial tsIdx cell order
fprintf('Initial recording depth estimates...\n')
% SubFunction:  Load a few thousand waveforms and use them to generate index sorted by depth of peak amplitude
[tsIdx, depthEst, depthOrder] = sortByDepth(tsIdx, kilopath, kiloInfo.chanMap);

% apply depth ordering to [sortId] (Kilosort generated)
spk.sortId = spk.sortId(depthOrder,:);
[~, nearestChan] = min(abs(depthEst-kiloInfo.chanMap.ycoords')');
% spk.id = nearestChan(:); % replacement/placeholder for former .id == (10*channel# + unit#);
% spk.id = cell2mat(spk.sortId(:,3)); % ...not super human-relevant, but 

spk.info = kiloInfo;
spk.sync = sync;

% strip out anything before sync trial start or more than a 2 minutes after last trial sync end.
% ...this could backfire, for some future extreme case.       --TBC Sept, 2014 --TBC, again in 2019


%%
if dosync
    % limit spike times to duration of experimental file (defined in .sync struct)
    tsSyncLimSamples = tsSyncLimSecs .* kiloInfo.sample_rate;
    tsIdx = cellfun( @(x) x(x>=tsSyncLimSamples(1) & x<=tsSyncLimSamples(end)), tsIdx, 'UniformOutput',0);
end


%% Get waveforms (or kilosort template approximations...or actual waveform means?)
%...matching rest of typical [spk] fields

% compile params struct for waveform extraction
clear pars
pars.nWf = 10000;                   % Number of waveforms per unit to pull out
pars.dataDir = kilopath;           % KiloSort/Phy output folder
% TODO: incorporate stereoprobe trodal waveform means & ci
pars.trodality = 1;% + contains(lower(spk.info.comment), 'stereo');
pars.chDepth = spk.info.chanMap.ycoords;

% % Most pars values have defaults best set in the loading function, but could be set manually
% pars.wfWin = [-40 41];             % Number of samples before and after spiketime to include in waveform
% pars.fileName = kiloInfo.dat_path;  % .dat file containing raw voltages (can be relative to .dataDir; e.g. '../../raw/file.dat')
% pars.dataType = kiloInfo.dtype;    % Data type of .dat file (this should be BP filtered)
% pars.nCh = kiloInfo.n_channels_dat;% Number of channels that were streamed to disk in .dat file

% gwfparams.spikeTimes = ceil(sp.st(sp.clu==155)*30000); % Vector of cluster spike times (in samples) same length as .spikeClusters
% gwfparams.spikeClusters = sp.clu(sp.clu==155);

[spk.wf, dupes, wfFigs] = getWaveforms_kilo(tsIdx, pars);

% Remove duplicate spikes (detected during waveform alignment & SNR calc)
tsIdx = cellfun(@(t,d) t(~d), tsIdx, dupes, 'uni',0);

% % Convert tsIdx from counts to seconds;   ts = {nUnits, 1}(nspikes,1))
% spk.ts = cellfun(@(x) x./kiloInfo.sample_rate, tsIdx, 'uni',0);   % ...NOTE: we still need tsIdx to load waveforms below

% Convert tsIdx from counts to seconds, AND sync seconds to PLDAPS frame time from using conversion function from .sync struct;
spk.ts = cellfun(@(x) sync.plx2ptb(x./kiloInfo.sample_rate), tsIdx, 'uni',0);   %.ts = {nUnits, 1}(nspikes,1))


%% Sync spikes with stimulus time
% Conversion functions are included in the .sync struct
%   -These aren't necessarily used if analysis all done relative to event syncs, but good to have around --TBC 2019


%% Compile output struct
% spike count for each unit
spk.n = cellfun(@length, spk.ts);
spk.snr = spk.wf.snr;


%% Depth estimate based on mean waveform amplitude
% peakCh already computed in getWaveforms_kilo.m
% ...no reordering here, just metrics
wfAmp = squeeze(range(spk.wf.mu));
% zero out small amp channels (mmany could pull center of mass)
% wfAmpThresh = 0.1*max(wfAmp);
wfAmpClipped = wfAmp;
wfAmpClipped(bsxfun(@lt, wfAmp, 0.2*max(wfAmp))) = eps;
wfDepth = (sum(bsxfun(@times, wfAmpClipped, spk.info.chanMap.ycoords))./sum(wfAmpClipped));

% [~, nearestChan] = min(abs(wfDepth(:)-kiloInfo.chanMap.ycoords')');
nUnits = length(spk.wf.peakCh);
spk.id = (1:nUnits)' + spk.wf.peakCh(:)./10^ceil(log10(nUnits)); % == unit#.nearestChan ...replacement/placeholder for former .id == (10*channel# + unit#);
spk.wf.depth = wfDepth(:);

% save wf figs as png only (qualitative diagnostic)
try
    if exist('wfFigs','var') %&& evalin('caller','saveout')
        for u = 1:length(wfFigs)
            figureFS(wfFigs{u},'portrait',[8 8]);
            drawnow
            [~,srcFile] = fileparts(sync.info.pdsSrc);
            set(gcf, 'name', sprintf('%s u%03d', srcFile,u));
            saveFigTriplet(0, sprintf('%s   %.2f  (=unit#.PeakChan#)', srcFile, spk.id(u)), [0,1,0,0,0], sprintf('wfFigs_kilosorted%s%s',filesep,srcFile));
        end
        % close all wf plot figs
        close([wfFigs{:}])
    end
end

%% --- end of main function commands ---


% % % % % % % % % %
%% Nested functions
% % % % % % % % % %


%% loadSortIds
    function [sortId, sortSrc] = loadSortIds(notThese, sortSrc)
        % Loads unit/cluster identity names as set by Kilosort or during manual curation in Phy
        % Returns:  [sortId] cell nUnits-by-3:  {sort#, sortIdString, kiloCluster#}
        %           [sortSrc] path to sortId file that was used
        %
        % No inputs necessary if everything is in default location in the kilopath.
        % [notThese] is a cellstr of sort id names that should be skipped
        %           def: {'noise'};
        % [sortSrc] is path to sort names...prob best left empty
        %           def: tries 'cluster_group.tsv', then 'cluster_groups.csv'
        
        % Parse inputs & defaults
        sortId = [];
        if nargin<1
            % sort ids to exclude
            notThese = {'noise'};
        end
        
        if nargin<2
            sortSrc = fullfile(kilopath, 'cluster_group.tsv');
            if ~exist(sortSrc, 'file')
                sortSrc = fullfile(kilopath, 'cluster_groups.csv');
                if ~exist(sortSrc, 'file')
                    sortSrc = [];
                    return
                end
            end
        elseif ~exist(sortSrc, 'file')
            warning('SortId file not found:\n\t%s\n', sortSrc);
            return
        end
        
        
        %% Read SortId source
        fid = fopen(sortSrc);
        C = textscan(fid, '%s%s');
        fclose(fid);
        % Format & scrap header line
        C = [lower(C{2}(2:end)), cellfun(@str2num, C{1}(2:end), 'uni',0)];
        
        
        %% Parse sortIds
        sortId = cell(size(C,1), 3);
        sortId(:,2:3) = C;
        
        % extract unique user-sepcified sort names
        uniSorts = unique([string(sortId(:,2)); "test01";"apples";"oranges"]);
        % exclude default names that have special handling already
        uniSorts = uniSorts(~contains(uniSorts, {'unsorted','good','mua','noise'}));
        
        % Assign sort names to standard values (klunky-ish for loop)
        for i = 1:size(sortId,1)
            switch sortId{i,2}
                case 'unsorted'
                    sortId{i,1} = 0;
                case 'good'
                    sortId{i,1} = 1;
                case 'mua'
                    sortId{i,1} = 2;
                case 'noise'
                    sortId{i,1} = -1;
                    
                otherwise
                    sortId{i,1} = 2 + find(uniSorts==sortId{i,2});
            end
        end
        
        sortId = sortId(~contains(sortId(:,2), notThese),:);
                
    end % loadSortIds


end %main function


% % % % % % % % %
%% Sub-functions
% % % % % % % % %


%% sortByDepth
function [tsIdx, depthEst, depthOrder] = sortByDepth(tsIdx, kilopath, chanMap)
%  Currently [tsIdx] & [depthEst] outputs are both presorted
%  	-Unsure if makes most sense to pass out:
%       index & depths(original)
%       index & depths(sorted)
%       tsIdx(sorted) & depths(sorted) & index
% 

% compile params struct for waveform extraction
%clear pars
pars.nWf = 2000;                   % Number of waveforms per unit to pull out
pars.dataDir = kilopath;           % KiloSort/Phy output folder
pars.trodality = 1;% + contains(lower(spk.info.comment), 'stereo');
% pars.chDepth = spk.info.chanMap.ycoords;

% % Most pars values have defaults best set in the loading function, but could be set manually
% pars.wfWin = [-40 41];             % Number of samples before and after spiketime to include in waveform
% pars.fileName = kiloInfo.dat_path;  % .dat file containing raw voltages (can be relative to .dataDir; e.g. '../../raw/file.dat')
% pars.dataType = kiloInfo.dtype;    % Data type of .dat file (this should be BP filtered)
% pars.nCh = kiloInfo.n_channels_dat;% Number of channels that were streamed to disk in .dat file

% gwfparams.spikeTimes = ceil(sp.st(sp.clu==155)*30000); % Vector of cluster spike times (in samples) same length as .spikeClusters
% gwfparams.spikeClusters = sp.clu(sp.clu==155);

wf = getWaveforms_kilo(tsIdx, pars, 'align',false, 'acg',false); % skip alignment, acg, & snr calc when just determining depth


%% Depth estimate based on mean waveform amplitude
wfAmp = squeeze(range(wf.mu));
[wfAmpMax, wfChPeak] = max(wfAmp);
[a, depthOrder] = sortrows([wfChPeak; wfAmpMax]',[1, -2]);

wfAmpClipped = wfAmp;
% % % % zero out small amp channels (mmany could pull center of mass)
wfAmpClipped(bsxfun(@lt, wfAmp, 0.2*max(wfAmp))) = eps;   %% wfAmpThresh = 0.1*max(wfAmp);
depthEst = (sum(bsxfun(@times, wfAmpClipped, chanMap.ycoords))./sum(wfAmpClipped));
% Order units from shallowest to deepest first, then largest to smallest amplitude second
%   (...SNR would be a better second dim, but not yet incorporated into Kilosort pipeline)
% depthSign = sign(diff(chanMap.ycoords([1,end])));
% [a, depthOrder] = sortrows([depthEst; wfAmpMax]',[depthSign, -2]);

% Sort outputs before returning
depthEst = a(:,1);
tsIdx = tsIdx(depthOrder);

end

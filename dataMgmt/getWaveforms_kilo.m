function wf = getWaveforms_kilo(tsIdx, varargin)
% function wf = getWaveforms_kilo(tsIdx, pars)
%
% Extracts individual spike waveforms from the raw datafile, for multiple
% clusters. Returns the waveforms and their means within clusters.
%
% Contributed by C. Schoonover and A. Fink
%
% % EXAMPLE INPUT
% tsIdx = cell of spike timestamps in samples;  {nUnits, 1}(nspikes,1))
% 
% pars.dataDir = '/path/to/data/';    % KiloSort/Phy output folder
% pars.fileName = 'data.dat';         % .dat file containing the raw 
% pars.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
% pars.nCh = 32;                      % Number of channels that were streamed to disk in .dat file
% pars.wfWin = [-40 41];              % Number of samples before and after spiketime to include in waveform
% pars.nWf = 2000;                    % Number of waveforms per unit to pull out
%     % pars.spikeTimes =    [2,3,5,7,8,9]; % Vector of cluster spike times (in samples) same length as .spikeClusters
%     % pars.spikeClusters = [1,2,1,1,1,2]; % Vector of cluster IDs (Phy nomenclature)   same length as .spikeTimes
%
% % OUTPUT
% wf.mu         % [nSamples, nCh, nUnits]       Average of all waveforms (per channel)
% wf.ci         % [nSamples, nCh, nCi, nUnits]  Conf intervals ( [.16,.84] prctiles)
%     % wf.unitIDs                               % [nClu,1]            List of cluster IDs; defines order used in all wf.* variables
%     % wf.spikeTimeKeeps                        % [nClu,nWf]          Which spike times were used for the waveforms
% wf.waveForms  % [nSamples, nCh, nWf, nUnits]            % [nClu,nWf,nCh,nSWf] Individual waveforms
%
% nUnits:   number of different units/clusters  == size(tsIdx,1) == numel(tsIdx) 
% nSamples: number of samples per waveform      == diff(pars.wfWin)
% 
% % USAGE
% wf = getWaveForms(pars);


%% Parse inputs & setup default parameters
pp = inputParser();
% Basic data parameters
pp.addParameter('dataDir', pwd, @isfolder)
pp.addParameter('fileName', []);
pp.addParameter('dataType', 'int16');
pp.addParameter('nCh', []); % Number of channels that were streamed to disk in .dat file
pp.addParameter('trodality', 1); % if >1, will reshape waveforms to n*wfSamples
pp.addParameter('chDepth', []); % depth of channel; input from spk.info.chanMap.ycoords (def 1:nCh)
pp.addParameter('sr', []);  % Raw data sampling rate
pp.addParameter('nWf', 1000);   % nWaveforms to load/average
pp.addParameter('wfWin', [-200, 800]*1e-6);   % waveform samples [before, after] timestamp (in seconds OR samples)
pp.addParameter('ci', [.16 .84]); % waveform ci levels (def +/-1sd)
pp.addParameter('all', false);   % return all waveforms loaded (i.e. not just mean & ci wfs)


%%
pp.parse(varargin{:});
%%
pars = pp.Results;

% Fill in the gaps
if isempty(pars.fileName)
    % load info from params.py
    kiloInfo = loadParamsPy( fullfile(pars.dataDir, 'params.py'));
    pars.fileName = kiloInfo.dat_path;  % .dat file containing raw voltages (can be relative to .dataDir; e.g. '../../raw/file.dat')
    pars.dataType = kiloInfo.dtype;     % Data type of .dat file (this should be BP filtered)
    pars.nCh = kiloInfo.n_channels_dat; % Number of channels that were streamed to disk in .dat file
    pars.sr = kiloInfo.sample_rate;
end
pars.datPath = fullfile(pars.dataDir, pars.fileName);

% Convert wf duration to samples
if abs(pars.wfWin(1))<1 || all(round(pars.wfWin)~=pars.wfWin)
    pars.wfWin = round(pars.wfWin*pars.sr);
end
if isempty(pars.chDepth), pars.chDepth = 1:pars.nCh; end

% create output struct
wf = struct('mu',[], 'ci',[], 'pars', pars);

if ~exist(pars.datPath, 'file')
    % Should have a standard .mat file of precompiled [wf] struct here to load
    % instead of always needing/reloading jumbo raw data
    fprintf(2, '~~~\tCould not find raw dat path to load waveforms from...continuing without them\n~~~\t%s\n', pars.datPath)
    return
end

%% Load .dat and KiloSort/Phy output
filenamestruct = dir(pars.datPath);
dataTypeNBytes = numel(typecast(cast(0, pars.dataType), 'uint8')); % determine number of bytes per sample
nSampTotal = filenamestruct.bytes/(pars.nCh*dataTypeNBytes);  % Number of samples per channel
datSize = [pars.nCh nSampTotal];
nSamples = length(pars.wfWin(1):pars.wfWin(end));

mmf = memmapfile(pars.datPath, 'Format', {pars.dataType, datSize, 'x'});
chMap = readNPY(fullfile(pars.dataDir, 'channel_map.npy'))+1;               % Order in which data was streamed to disk; must be 1-indexed for Matlab
nChInMap = numel(chMap);
nUnits = numel(tsIdx);

%%
% Read spike time-centered waveforms
%     unitIDs         = unique(pars.spikeClusters);
%    nUnits        = size(unitIDs, 1);
% spikeTimeKeeps  = nan(nUnits, pars.nWf);

%%
clear spk
% Slice variables for parfor
[allWfs, tmpMu, tmpCi] = deal(cell(nUnits,1));
nWf = pars.nWf;
nCh = pars.nCh;
wfWin = pars.wfWin(1):pars.wfWin(end);
ciLvls = pars.ci;
% Only compile ALL waveforms if requested (big & slow!)
getAllWfs = pars.all;

tic
fprintf('\nTotal units:\t[%s]\nLoading wfs:\t[', repmat('.',[nUnits,1]))
parfor u = 1:nUnits
    fprintf('.')
    tsU = tsIdx{u};
    theseWf = min([nWf, length(tsU)]);
    % t0 index for each waveform
    wf0 = tsU(randperm(length(tsU), theseWf));    % spikeTimesRP = curSpikeTimes(randperm(curUnitnSpikes, nWf);
    %...thus, all wf timepoints ==
    wfi = (wf0 + wfWin)'; % wf index of [nSamples, nWf]       %   spikeTimeKeeps(u,1:nWf) = sort(spikeTimesRP(1:nWf));
    % Put indices in wfs variable (for efficiency)
    tmpWf = mmf.Data.x(1:nCh, wfi);     %#ok<PFBNS>
    tmpWf = permute( reshape(tmpWf(chMap,:), [nCh, size(wfi)]), [2,1,3]);
    tmpMu{u} = nanmean(tmpWf, 3);
    tmpCi{u} = quantile(tmpWf, ciLvls, 3);
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
% wf.mu          = nan(nSamples, pars.nCh, nUnits);
% wf.ci          = nan(nSamples, pars.nCh, numel(pars.ci), nUnits);

% Only deal with all wfs if we have to...
if getAllWfs
    wf.all    = nan(nSamples, pars.nCh, pars.nWf, nUnits);
    for u = 1:nUnits
        wf.all(:,:,1:size(allWfs{u},3),u) = allWfs{u};
    end
end

% % Compute summary stats
% wfAmp = squeeze(range(wf.mu), 1);
% [wf.stat.peakAmp, wf.stat.peakCh] = max(wfAmp);
% % zero out small amp channels (mmany could pull center of mass)
% wfAmpClipped = wfAmp;
% wfAmpClipped(bsxfun(@lt, wfAmp, 0.1*max(wfAmp))) = 0;


fprintf('\tdone. (%3.1f sec)\n', toc)


%%

% Complete spk output struct

end

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

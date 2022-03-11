function plxInfo = getPlxInfo(plxFilename, plotEvents)
%
% Retrieve all info about PLX/PL2 file [plxFilename]% 
% -- A helper function to simplify extraction of all outputs from the two
% info functions in the Plexon SDK (plx_info.m & plx_information.m)
% 
% Returns plxInfo struct with fields:
%     [path]
%     [ver]
%     [freq]
%     [comment]
%     [trodality]
%     [wfSize]
%     [preThresh]
%     [spikePeakV]
%     [spikeBitRes]
%     [slowPeakV]
%     [slowBitRes]
%     [duration]
%     [dateStr]
%     [adNames]
%     [adNumCountFreqGains]
%     [tscounts]        (plx_info.m)
%     [wfcounts]        (plx_info.m)
%     [evcounts]        (plx_info.m)
%     [contcounts]      (plx_info.m)
%     [strobes]
%     [rstart]
%     [rstop]
%     [plx2clock]
% 
% 2018-12-02  TBC  Wrote it.
% 2021-01-08  TBC  Added basic event loading (for determining where experiment files were run w/in longer plx files)
%                  .plx2clock converts plx event times to [approx] real clock time (for comparison to PLDAPS file times)
% 2022-01-26  TBC  Kilosort workaround: Allow plxInfo struct to be loaded from saved rawInfo.mat file
%                  - Workaround for manual merging pl2 files 


%% Defaults
if nargin<2 || isempty(plotEvents)
    plotEvents = 0;
end


%% General PLX info struct
tic
fn = {'path', 'ver', 'freq', 'comment', 'trodality', 'wfSize', 'preThresh', 'spikePeakV', 'spikeBitRes', 'slowPeakV', 'slowBitRes', 'duration', 'dateStr'};
nfo = cell(1,length(fn));

try
    % load PLX file info directly
    [nfo{:}] = plx_information(plxFilename);
catch
    % if fail, attempt to load plxInfo struct
    load(plxFilename, 'plxInfo', '-mat');
    return
end

plxInfo =  cell2struct(nfo, fn, 2);

% Additional pl2 file details
[ ~, isPl2 ] = internalPL2ResolveFilenamePlx( plxFilename );
if isPl2
    pl2 = PL2GetFileIndex(plxFilename);
    plxInfo.CreatorSoftwareName = pl2.CreatorSoftwareName;
    plxInfo.CreatorSoftwareVersion = pl2.CreatorSoftwareVersion;
    if strcmpi(pl2.CreatorSoftwareName(1:2),'SC')
        % "Sort Client" ...retro control software for Plexon MAP system
        % - additional checks & 
        plxInfo.isMAP = true;
    else
        plxInfo.isMAP = false;
    end
    
end


% AD channel info (continuous spike & LFP)
[adNames, adRawNum, adCount, adFreq, adGains] = deal([]);
try % prevent crash if no AD channels
    [n, adNames]    = plx_adchan_names(plxFilename);
        adNames     = cellstr(adNames);
    [n, adRawNum]   = plx_ad_chanmap(plxFilename);
    [n, adCount]    = plx_adchan_samplecounts(plxFilename);
    [n, adFreq]     = plx_adchan_freqs(plxFilename);
    [n, adGains]    = plx_adchan_gains(plxFilename);
catch
    fprintf(2, '\tProblem while reading analog data info from plx file.\n\tLikely just buggy SDK functions. Carry on, but be warned...\n')
end


%% Discard empty channels and low sampling rate continuous channels (LFP)
ii = adCount>0 & adFreq==max(unique(adFreq));

plxInfo.adNames = adNames(ii);
plxInfo.adNumCountFreqGains = [adRawNum(ii)', adCount(ii), adFreq(ii), adGains(ii)];

% Channel count (explicitly)
if isPl2
    pl2 = PL2GetFileIndex(plxFilename);
    plxInfo.nChannels = pl2.NumberOfRecordedSpikeChannels;
else
    plxInfo.nChannels = size(plxInfo.adNames,1);
end

% slightly more nitty-gritty counts from plx_info.m
[plxInfo.tscounts, plxInfo.wfcounts, plxInfo.evcounts, plxInfo.contcounts] = plx_info(plxFilename);


%% get PLX events
% try
    plxInfo = getPlxEvents(plxInfo, plotEvents);
% end


%% Done.
% report to command window
% fprintf('\t%2.2f sec to query info on plx file:\t%s\n', toc, plxFilename);


end %main function

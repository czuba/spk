function plxInfo = getPlxInfo(plxFilename)
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
% 
% 2018-12-02  TBC  Wrote it.

% General PLX info struct
tic
fn = {'path', 'ver', 'freq', 'comment', 'trodality', 'wfSize', 'preThresh', 'spikePeakV', 'spikeBitRes', 'slowPeakV', 'slowBitRes', 'duration', 'dateStr'};
nfo = cell(1,length(fn));

[nfo{:}] = plx_information(plxFilename);
plxInfo =  cell2struct(nfo, fn, 2);

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
ii = adCount>0; % discard empty channels
plxInfo.adNames = adNames(ii);
plxInfo.adNumCountFreqGains = [adRawNum(ii)', adCount(ii), adFreq(ii), adGains(ii)];

% slightly more nitty-gritty counts from plx_info.m
[plxInfo.tscounts, plxInfo.wfcounts, plxInfo.evcounts, plxInfo.contcounts] = plx_info(plxFilename);


fprintf('\t%2.2f sec to query info on plx file:\t%s\n', toc, plxFilename);
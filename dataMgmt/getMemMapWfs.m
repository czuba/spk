function wfs = getMemMapWfs(mmf, ts, wfWin, chans)
% function wfs = getMemMapWfs(mmf, ts, wfWin, chans)
% 
% Gets waveform snippets from memory mapped binary file[mmf]
% INPUTS:
% [mmf]     Handle to memory mapped file; mapped as mmf.Data.x(channels, samples)
%           % mmf   = memmapfile(datPath, 'Format', {'int16', [nchan, totSamples], 'x'});
% [ts]      Spike timestamps
% [wfWin]   Window of data samples to be read, relative to each timestamp[ts] 
%           % wfWin = [-20, 40];
% [chans]   Channel indices 
% 
% OUTPUTS:
% [wfs]     Waveforms with dimensions:  [wfSample, channel, count]
% 
% 2020-08-xx TBC Wrote it.
% 

% expand timestamps to wf sample indices
if numel(wfWin)==2
    wfWin = wfWin(1):wfWin(end);
elseif numel(wfWin)==1
    wfWin = -wfWin:wfWin;
end

nch = numel(chans);
% Put sample indices in wfs variable (for efficiency?)
wfi = (ts + wfWin)'; % wf index of [nSamples, nWf]
tmpWf = mmf.Data.x(chans, wfi);
% reshape as [wfSample, channel, count]
wfs = permute( reshape(tmpWf, [nch, size(wfi)]), [2,1,3]);

end %main function

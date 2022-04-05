function [wf2, lag] = alignWfs_kilo(wf, wEnd, maxShift, nPass)
% function [wf2, lag] = alignWfs_kilo(wf, wEnd, maxShift, nPass)
% 
% Aligns waveforms by estimating the best fitting delay/shift (crosscorrelation via finddelay.m)
%   [wf]        = Matrix of waveforms: [count, nsamp]  (matching returned by plx_waves)
%   [wEnd]      = Meant to clip wf length
%                   ~!~ NOT implimented ~!~  Defaults to all points in wf: size(wf,2)   TBC 2014-12-01
%                   --minimize noise/var from wf length exceeding spike width (e.g. 1200µs) and clipped wf tails
%                   --ideally would really be range of wf points to use, but implementation is shifty b/c of device
%                     dependent sampling rate...
%   [maxShift]  = Maximum waveform shift (in wf samples [unfortunately] NOT ms)
%                   --Default == 10
%   [nPass]     = Iterative alignment passes through data (def = 1)
% 
% 2014-09-XX TBC  Wrote it
% 2014-12-01 TBC  Commented
% 2020-08-18 TBC  Updating for kilosort2 waveform samples

[nwf, nsamp] = size(wf);

if ~exist('wEnd','var') || isempty(wEnd)
    wEnd = nsamp;
end
if ~exist('maxShift','var') || isempty(maxShift)
    maxShift = 40;
end
if ~exist('nPass','var') || isempty(nPass)
    nPass = 1;
end

% Must do coarse alignment before computing mean
%   - meant to address clusters with temporal offsets due to triggering/template fit issues
%   - ...this is soo hairy/ugly.
nmu = min([nwf, 2000]);     % use up-to 2,000 wfs for mean
minPt = floor(nsamp/3);     % align troughs 1/3 into wf sample range
wfi = randperm(nwf, nmu);   % randomly sample from input wfs
[~, minAmp] = min(wf(wfi,1:2*minPt), [],2); % find troughs within first 2/3 of waveform

minAmp = minAmp-minPt;  % make relative to desired trough alignment
wfmu = nan(nmu, nsamp+2*minPt);
for i = 1:length(wfi)
    ii = (1:nsamp)+minPt-minAmp(i); % trough-aligned insertion point
    wfmu(i, ii) = wf(wfi(i),:);
end
wfmu = mean(wfmu(:, minPt+1:end-minPt),'omitnan');
wfmu(isnan(wfmu)) = 0;

% find align shifts
lag = -finddelay(wf', repmat(wfmu,[nwf,1])', maxShift)';

% if mean(abs(lag))>5
%     keyboard
% end

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

% Repeat as requested
if nPass>1
    [wf2, lag2] = alignWfs_kilo(wf2, [], maxShift, nPass-1);
    % accumulate offsets
    lag = lag+lag2;
end

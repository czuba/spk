function snr = calcWaveformSNR(mmf, wf0, wfWin, chan, makePlots)
% function snr = calcWaveformSNR(mmf, wf0, wfWin, chan, makePlots)
% 
% Uses memmapped file handle [mmf] to read waveforms of size [wfWin] from each timestamp [wf0]
% Computes spike waveform [snr] from gaussian approx. of projection across all input [chan]
% 
% NOTE: important to align waveform timestamps prior to snr calc,
% else timeshifted units likely appear artificially low-SNR
% 
% [makePlots]   default==true
% 
% see also: getMemMapWfs
% 

if nargin<5 || isempty(makePlots)
    makePlots = 1;
end

% wf0 are spike timestamps
wfs = getMemMapWfs(mmf, wf0, wfWin, chan);


% Compute mean waveform
% NOTE: important to align waveform timestamps prior to snr calc,
% else timeshifted units will appear artificially low-SNR
mnWF = mean(wfs, 3); 

% for projections, reshape and just multiply
mnWFlin = reshape(mnWF,1,[]);
pOwn = mnWFlin*double(reshape(wfs, numel(mnWF), []));

if makePlots
    % plot mean waveform
    H = figure; 
    subplot(3,2,1);
    if min(size(mnWF))<3
        plot(mnWF)
        ylabel('meas V');
    else
        imagesc(mnWF');
        ylabel('channel');
    end
    fmtPlot('small');
    title('mean waveform');
    %xlabel('samples'); 
    
    subplot(3,2,2);
    if min(size(mnWF))<3
        plot(bsxfun(@rdivide, mnWF, max(abs(mnWF),[],2)+3))
    else
        imagesc(bsxfun(@rdivide, mnWF, max(abs(mnWF),[],2)+3)');
        ylabel('channel');
    end
    fmtPlot('small');
    title('rescaled');
%     drawnow;
end

%% 2. read many other random points and get their projections onto the mean

null0 = randi([min(wf0),max(wf0)], length(wf0),1);

wfNull = getMemMapWfs(mmf, null0, wfWin, chan);
pOther = mnWFlin*double(reshape(wfNull, numel(mnWF), []));
% % %     % now pick random times to look at - but NOT random times around existing
% % %     % spikes. 
% % %     sampsTaken = ceil(spikeTimes*datPars.Fs)+(datPars.wfWin(1):datPars.wfWin(end));
% % %     allSamps = 1:nSamp;
% % %     sampsAvail = allSamps(~ismember(allSamps, sampsTaken(:))); 
% % %     wfSamps = sampsAvail(randi(numel(sampsAvail), [datPars.nSpikesToUse 1])); 
% % %     wfs = zeros(numel(datPars.chanMap), wfNSamps, numel(wfSamps));
% % %     for q = 1:numel(wfSamps)
% % %         tmpWf = mmf.Data.x(:,wfSamps(q):wfSamps(q)+wfNSamps-1);
% % %         wfs(:,:,q) = tmpWf(datPars.chanMap,:);
% % %     end
% % %     pOther = mnWFlin*reshape(wfs, numel(mnWF), []);

%% 3. plot histogram/density, compute SNR by gaussian approximation

snr = (mean(pOwn)-mean(pOther))./std(pOther);

if makePlots
    sp = subplot(3,2,4);
    bins = linspace(min([pOwn pOther]), max([pOwn pOther]), 101);
    [n,x] = hist(pOwn, bins); 
    plot(sp, x,n); hold on; 
    [n,x] = hist(pOther, bins); 
    plot(sp,x,n);
    title(sp,sprintf('snr = %.2f', snr))
    xlabel(sp,'projection onto mean');
    fmtPlot('small');
    lh = legend({'spikes', 'other'}, 'fontsize',6, 'box','off');
    drawnow
end


end %main function

function H = addWfInset(AH, wf, u)

% create new inset axes
% % % try
% % %     axpos = plotboxpos(AH); % Unpack Matlab's convoluted position params
% % %     AHi = axes('units','normalized','position', [axpos(1)+1.1*axpos(3), axpos(2)-0.5*axpos(4), 0.8*axpos(3), 0.8*axpos(4)]);
% % % 
% % % catch
    axpos = get(AH, 'position');
    AHi = axes('position', [axpos(1)+0.7*axpos(3), axpos(2)-0.25*axpos(4), 0.4*axpos(3), 0.2*axpos(4)]);
% % % % AHi = axes('position', [axpos(1)+0.9*axpos(3), axpos(2), 0.4*axpos(3), 0.4*axpos(4)]);
% % % end
% slight tweak ax color so easily selectable in post-production
axcol = .07*[1 1 1];    set(AHi, 'xcolor',axcol, 'ycolor',axcol)

if ndims(wf.mu)==3
    % Only plot the wf channel with the highest amplitude
    wfVals = wf.mu(:,:,u);
    [~,ii] = max(range(wfVals));
    wfVals = wfVals(:,ii);
    wfValsCi = squeeze(wf.ci(:,ii,:,u));
else
    wfVals = wf.mu(:,u);
    wfValsCi = wf.ci(:,:,u);
end

% WF conf interval in red(ish)
plot(AHi, wfValsCi, 'color',[1,axcol(2:3)], 'linewidth',.5);
hold on,
% WF mean in black(ish)
plot(AHi, wfVals, 'color',axcol, 'linewidth',.5);

if range(wf.mu(:))<1
    % wf vals in mV
    yl = max(0.8*range(wf.mu(:)), .1);     % yl = .1; %yl = .04;
else
    % else...wf vals in arb. a/d units (3000 is good)
    yl = min(0.8*range(wf.mu, 'all'), 800); %3000;
end
yt = yl/2;
set(AHi, 'fontsize',6);
set(AHi, 'color','none', 'yAxisLocation','right');
set(AHi, 'linewidth',0.5, 'plotboxaspectRatio',[3,2,1],...
    'yticklabel',[], 'ytick',yt*[-1,0,1], 'ylim',yl*[-.5,.5],...
    'xlim',[1,size(wf.mu,1)], 'xticklabel',[]);
box off

if ndims(wf.mu)==3
    % Kilosort output with nCh waveforms per unit
    set(AHi, 'ytick', wfVals(end), 'yticklabel',ii);
    %     % don't have access to this in .wf struct.   TODO: fix after fixing [main] Kilosort2 issues
    %     xlabel(dv.uprb.sortId{u,2})
end
uistack(AHi, 'top')
if nargout==0
    return
else
    H = AHi;
end

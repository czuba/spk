function H = addWfInset(AH, wf, u, relShift, sc)
%

if nargin<5 || isempty(sc)
    % no scaling by default
    sc = 1;
elseif isa(sc,'struct')
    % standard "dv.uprb.info" struct
    if isfield(sc,'dat_path') && isfield(sc,'rawInfo')
        % pull conversion fxn from inside [messy] hierarchy
        sc = sc.rawInfo.rawInfo.raw2v;
    end
end
% Continuous raw (int16) data to mv
if isa(sc,'function_handle')
    sc = sc(1);
end

if nargin<4 || isempty(relShift)
    relShift = [0 0];
end

% create new inset axes
if ~isscalar(AH) && numel(AH)==3
    % use as inputs to subplot
    AHi = subplot(AH(1), AH(2), AH(3));
    isInset = false;
    
elseif ishandle(AH)
    if contains(class(AH),'matlab.graphics.axis') %isa(AH,'matlab.graphics.axis.Axes')
        % "if contains(..." works for polarAxes too
        % get axes position
        axes(AH);
        axpos = get(AH, 'position');
        % bump over & up a bit more for polar plots
        if contains(class(AH), 'Polar')
            relShift = relShift + [.1, -.35];
        end
        % tediously define location for axes inset
        AHi = axes('position', [axpos(1)+(0.7+relShift(1))*axpos(3), axpos(2)-(0.25+relShift(2))*axpos(4), 0.4*axpos(3), 0.2*axpos(4)]);
        % % % % AHi = axes('position', [axpos(1)+0.9*axpos(3), axpos(2), 0.4*axpos(3), 0.4*axpos(4)]);
        isInset = true;
        
    elseif isa(AH,'matlab.ui.Figure')
        % if figure handle
        figure(AH);
        AHi = findobj(AH,'type','axes');
        if ~isempty(AHi)
            % add inset to most recent axis
            AHi = AHi(1);
            isInset = true;
        else
            % or create first axes in figure
            AHi = axes(AH);
            isInset = false;
        end
    end
else
    return % silently move on
end

% slight tweak ax color so easily selectable in post-production
axcol = .07*[1 1 1];    set(AHi, 'xcolor',axcol, 'ycolor',axcol)

if ndims(wf.mu)==3
    % Only plot the wf channel with the highest amplitude
    wfVals = wf.mu(:,:,u);
    if isfield(wf, 'peakCh')
        ii = wf.peakCh(u);
    else
        [~,ii] = max(range(wfVals));
    end
    wfVals = wfVals(:,ii) .* sc;
    wfValsCi = squeeze(wf.ci(:,ii,:,u)) .* sc;
else
    wfVals = wf.mu(:,u) .* sc;
    wfValsCi = wf.ci(:,:,u) .* sc;
end

% WF conf interval in red(ish)
plot(AHi, wfValsCi, 'color',[1,axcol(2:3)], 'linewidth',.5);
hold on,
% WF mean in black(ish)
plot(AHi, wfVals, 'color',axcol, 'linewidth',.5);

wfrng = range(wf.mu(:) .* sc);
if wfrng<1
    % wf vals in mV
    yl = max(0.8*wfrng, .1);     % yl = .1; %yl = .04;
    % adjust yticks, round to 1 sigfigs
    yt = round((yl/2)*[-1,0,1], 1,'significant');
else
    % else...wf vals in arb. a/d units (3000 is good)
    yl = min(0.8*wfrng, 800); %3000;
    % adjust yticks, round to 2 sigfigs
    yt = round((yl/2)*[-1,0,1], 2,'significant');
end

% Axes formatting
box off
set(AHi, 'color','none', 'Clipping','off');
if isInset
    set(AHi, 'fontsize',6, 'yAxisLocation','right');
else
    fmtPlot('small')
    title('Peak channel waveform')
    set(AHi, 'fontsize',10);
    % if not inset, add legend
    %lh = legend({'ci','mean'}, 'fontsize',6, 'box','off');
end
set(AHi, 'linewidth',0.5, 'plotboxaspectRatio',[3,2,1],...
    'ytick',yt, 'ylim',yl*[-.55,.55],...
    'xlim',[1,size(wf.mu,1)], 'xticklabel',[], 'ticklength',[.03 .015]);

if ndims(wf.mu)==3
    % Kilosort output with nCh waveforms per unit
    ytl = vec2tick(yt);
    ytl{2} = sprintf(' ch%d',ii);
    if isInset
        ytl{1} = ''; ytl{3} = '';
    end
    set(AHi, 'ytick',yt, 'yticklabel',ytl);
end
if isfield(wf,'snr')
    xlabel(sprintf('snr = %.2f',wf.snr(u)));
end

uistack(AHi, 'top')
if nargout==0
    return
else
    H = AHi;
end

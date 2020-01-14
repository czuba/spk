function [dirpref, dirfit] = plotUnitDirTune(dv, u, sp)
% function [dirpref, dirfit] = plotUnitDirTune(dv, u, sp)
% 
% Drop-in function to add a direction tuning [sub]plot to a figure.
% (...presumably to relate tuning to a bunch of other subplots depicting other tuning/selectivity)
% 
% INPUTS:
%   [dv]    "dv" struct from dirTune_ file (typ generated w/ glDraw branch)
%   [u]     unit to plot (index, NOT 'unitID', e.g. ch12,unit1==unitId 121)
%   [sp]    handle to desired subplot; will plot into current axes if nothing input
% 
% OUTPUTS:
%   [dirpref]   vector average tuning [pref, bias]
%   [dirfit]    struct with fitted fxn handle & parameters (if available).
%               Eval with:  xs = -180:180; ys = dirfit.fxn(xs, dirfit.fit); plot(xs,ys)
% 
% 2019-11-20  TBC  Wrote it based on tuneXzAwake [...probably]

% extract local vars (copy-pasta from initial fitting fxn)
fileName = dv.info.pdsName;
stimType = dv.info.stimType;
dirs = dv.tune.dirs;



%% Figure & Plotting params
% fprintf('Plotting direction tuning')

msz = 8;
lw = 1.52;%1.8;

doPolar = false; % quite possible this plotting format branch is broken/old

cols = get(0,'DefaultAxesColorOrder');
refcol = .91*[1 1 1];
xt = [-90:90:450];
xtl = compose('%0.5g', xt);
xtl(1:2:end) = {''};


% Stimulus specifics
switch stimType{1}
    case 'gabors'
        dirs = wrapTo360(dv.tune.dirs); dirs(dirs==360) = 0;
        % Only use 2d directions (...ignore twd & away conds for now)
        dd = dirs(:,1)==dirs(:,2);
        dirs2d = dirs(dd,1);
        
    case 'dotBall'
        dd = dirs(:,2)==0;
        dirs2d = dirs(dd,1);
        
end

% dir range for tuning curve plots
if doPolar
    dirRng = 0:1:360;
    
else
    dirRng = -10:1:370;
    [dirsXax, dirsII] = sort(dirs2d, 'ascend');
    
    xl = [min(dirRng), max(dirRng)] + [-15,15]; %    xl = [min(dirRng), max(dirRng)+15];
end


%% Do the plotting

% 'plotBoxAspectRatio'
pbaxr = [5,3,1];

if ~exist('sp','var') || isempty(sp)
    sp = axes('plotboxaspectratio',pbaxr);
    %     % subplot spacing ([a.u.])
    %     spMargin = [.02 .02];
    %     sp = subplot_tight(spy ,spx, u+spbase, spMargin, 'plotboxaspectratio',pbaxr);
else
    cla(sp); 
    set(sp, 'plotboxaspectratio',pbaxr);
end

spPos = get(sp,'position');
set(sp, 'colororderindex', 1)

% janky fix for speed dimension
if ndims(dv.tune.trMu)==2
    jjj = 1;
else
    jjj = 1:size(dv.tune.trMu,2);
end

for ui = jjj
    set(sp, 'colororderindex', ui)
    trTmp = dv.tune.trMu(dd,ui, u);
    [db, dp, ob, op] = orivecfit(dirs2d, trTmp);
    
    % tuning plot limits
    yl = [-0.1, 1.1] *round2( max(trTmp)+2.5, 5);
    if yl(2)<30, yl(2)=30; end  % set ylim >= 30 spk/s
    rticks = unique([0, round2( max(trTmp(:))*[0.5,1], 5)]);
    
    % Plot on polar axes, wrapping endpoint
    if doPolar
        plh = polarplot(d2r(dirs2d([end,1:end])), trTmp([end,1:end]),'.-');
        sp = gca; set(sp,'position',spPos); hold on;  % srsly, wtf!?...quirky Matlab complication makes axes and polaraxes non-interchangeable
        curveColor = get(plh, 'color');
        
        set(sp, 'ThetaTick',[0:45:360], 'RTick', rticks)
        
    else
        % plot data
        plh = plot(sp, dirsXax, trTmp(dirsII), '-');
        hold on, box off
        curveColor = get(plh, 'color');
        
        plh(end+1) = plot(sp, dirsXax, trTmp(dirsII), '.');
    end
    
    set(plh, 'color', curveColor, 'linewidth',lw, 'markersize',msz);
end


xlabl = sprintf('ch.%d:\nO(%2.0f, %1.2f)    D(%2.0f, %1.2f)', dv.uprb.id(u), [op,ob], [dp,db]);

set(sp, 'linewidth',0.5)

if doPolar
    set(sp, 'ThetaTick',thTic, 'ThetaTickLabel', thTicl, 'RTick', rticks, 'rlim',yl);
else
    set(sp, 'plotBoxAspectRatio',pbaxr, 'XTick',xt,'XTickLabels',xtl,  'YTick', [rticks]);
    ylim(sp, yl); % [-15,375])
    xlim(sp, xl);   %[-15,375]);
end


titext = {fileName, sprintf('DirTune spd: %s deg/s', mat2str(dv.tune.spd))};

ht = title(titext, 'fontsize',10, 'fontweight','normal', 'interpreter','none');
xlabel(xlabl, 'fontsize',10, 'fontweight','normal')

% % draw stim loc on rf plot
% try % prevent error crash if no associated rf plot
%     stimPos = [commonParams.pos(1:2), commonParams.sz];
%     ha = findobj(figure(rfFigH), 'type','axes','tag',sprintf('ch%drf',u));
%     if ~isempty(ha)
%         set(ha, 'NextPlot','add');  %axes(ha); hold on;
%         hacl = get(ha,'clim');
%         
%         % Draw contours for gabor contrast envelope at [0.2, 0.5, 0.95] intensity
%         gw = commonParams.sz; % stim fwhm    tstStim(3)/2.355;
%         gww = gw*1.2; gx = linspace(-gww,gww,100); [gx,gy] = meshgrid(gx, gx); gg = zeros(size(gx)); gg(:) = hypot(gx(:), gy(:));   % gauss2D_R(gx, gy, gw, gw, 0, 1);
%         contour(ha, gx+stimPos(1), gy+stimPos(2), gg, gw*[1 1], 'color',curveColor(end,:), 'linewidth',1);
%         
%         set(ha, 'clim',hacl);
%         drawnow nocallbacks
%     end
% end

% Add waveform inset (to subplot[sp])
% ...this should really be its own compartmentalized fxn.   2019: IT IS!!
if ~isempty(dv.uprb.wf.mu)
    addWfInset(sp, dv.uprb.wf, u);
end


%% output tuning params
dirpref = [dp, db];
if isfield(dv.tune,'fxn')
    % fit params
    dirfit.fxn = dv.tune.fxn;
    dirfit.fit = dv.tune.fit(:,u);
    
else
    dirfit  = struct('fxn',[],'fit',[]);
end

% fprintf('\n')

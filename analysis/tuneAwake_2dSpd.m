function [dv] = tuneAwake_2dSpd(filebase, unitArgs, saveout, rfFigH)

%
% Modified version to drop out 3d motion conditions that were not presented in stereo
%


figDir = fullfile(pwd,'figs'); %
evalin('caller', sprintf('figDir = ''%s''', figDir))

%% Default inputs

doPolar = 0;

% file name
if nargin<1 || isempty(filebase)
    filebase = []%
end

% spikes to load/use
defarg('unitArgs', [666, 0]) % Default:  666=load all spikes, collapse all sorted & unsorted

if length(unitArgs)<2
    unitArgs = [unitArgs, 0];   % 2nd arg wf loading: 0=no[default], 1=yes & return ch stats,2=yes & return all wfs)
end

% save figs
defarg('saveout',nan)

if saveout>0 || ishandle(saveout)
    rfFigH = saveout;
elseif nargin<4 || isempty(rfFigH)
    rfFigH = [];
else  % limit to one linked rf figure
    rfFigH = rfFigH(1);
end

if ~exist('spkSrc','var')
    spkSrc = cellstr('uprb');
end


%% Load the files
if isstruct(filebase)
    % first input was existing dv struct, not filebase string
    dv = filebase;
else
    dv = initDaily(filebase, pwd, spkSrc, unitArgs);
    
    
    % Make sure this is data from the proper stimulus file type
    if isfield(dv.pds.baseParams.session, 'caller') && ~contains(lower(dv.pds.baseParams.session.caller.name), 'dirtune')
        warning('PDS file generated from:  %s\nThis doesn''t appear to be the correct file type for computing Direction Tuning.', dv.pds.baseParams.session.caller.name)
        keyboard
    end
    
end

fileName = dv.pds.baseParams.session.file(1:end-4);


%% Compute rate from condMatrix

% calc spike rate for each presentation
rate = calcRate_condMatrix(dv.uprb, dv.pds);

nunits = size(rate.count,2);


%% Parse condMatrix conditions & collapse across nessary dimensions
% Dependent on structure of ACTUAL MATRIX MODULE    (...in this case:  glDraw.doDirTune_gab.m)
baseModule = dv.pds.condMatrix.modNames.currentStim{1};

% Type of Tuning stimulus used (e.g. 'gabors' or 'dotBall')
stimType = {baseModule, dv.pds.baseParams.session.caller.name};

% --!!-- DANGER::  Assumes these params are FIXED throughout this file
switch stimType{1}
    case 'gabors'
        commonParams.pos = dv.pds.baseParams.(baseModule).stimCtr;
        commonParams.pos(2) = -commonParams.pos(2); % flip pos y-dimension when rendering texture stimuli
        %*** (...fix this in stim code to make universally consistent: negative below horizontal meridian)
        commonParams.sz = dv.pds.baseParams.(baseModule).gabSz;
        commonParams.sf = dv.pds.baseParams.(baseModule).gabSf;

        nDirs = rate.condDims(1);
        nSpd = rate.condDims(2);
        matSz = [nDirs, nSpd];
        
        % CONDITIONS:  [ LeftEyeDir, RightEyeDir, speedDegPerSec, posXpx, posYpx ]
        condPars = {'LeftEyeDir', 'RightEyeDir', 'speedDegPerSec', 'posXdeg', 'posYdeg'};
        condMat = cellfun(@(x) [x.dir, x.gabTf ./ commonParams.sf, x.stimPos.*[1,-1] + commonParams.pos(1:2)], dv.pds.condMatrix.conditions, 'uni',0);
        condSet = cell2mat(condMat(:));
        
        % Stim conditions for ALL strobes (in order of presentation)
        strobeConds = condSet(rate.stimStrobes(:,2),:);
        
        [diri, axi, oi] = ind2sub(rate.condDims, rate.stimStrobes(:,2));
        [xyi] = sub2ind(matSz, diri, axi); % ...this is only used as a counter for trial reps at each dir
        
        % collapse across ori/direction
        [dirs] = unique(condSet(:,1:2), 'rows','stable');
        [spd] =  unique(condSet(:,3),'stable');

        
    case 'dotBall'
        commonParams.pos = dv.pds.baseParams.(baseModule).stimCtr;
        commonParams.gridSz = dv.pds.baseParams.(baseModule).gridSz;
        commonParams.sz = dv.pds.baseParams.(baseModule).ballSz;    %gabSz;
        commonParams.spd = dv.pds.baseParams.(baseModule).dotSpd;
        nDirs = rate.condDims(1);
        nAxes = rate.condDims(2);
        nSpd = length(commonParams.spd);   %
        
        % Only frontoparallel directions for starters
        matSz = [nDirs, nAxes];
        
        % CONDITIONS:  [ xyDir, motPlane, speedDegPerSec, posXpx, posYpx ]
        condMat = cellfun(@(x) [x.dir, x.stimPos(1:2).*commonParams.gridSz(1:2) + commonParams.pos(1:2)], dv.pds.condMatrix.conditions, 'uni',0);
        condSet = cell2mat(condMat(:));
        
        % Stim conditions for ALL strobes (in order of presentation)
        strobeConds = condSet(rate.stimStrobes(:,2),:);
        
        % Convert strobe values into x & y indices, while collapsing across any
        % additional condMatrix dimensions (e.g. ori, direction, disparity, speed etc)
        
        % -- final [oi] output from ind2sub captures all [o]ther dimensions
        [diri, axi, oi] = ind2sub(rate.condDims, rate.stimStrobes(:,2));
        [xyi] = sub2ind(matSz, diri, axi); % ...this doesn't apply right now
        is2d = axi==1;
        
        % Set of unique directions
        dirs = condSet(unique(xyi), 1:2);
        ds = dirs(1:rate.condDims(1), 1);       % directions
        as = dirs(1:rate.condDims(1):end, 2);   % axes (??)
        
        
end


% trial count in each condition
ntr = reshape(hist(xyi, prod(matSz)), matSz);
% ntr = reshape(hist(diri(is2d), length(dirs)), matSz);

[ct tr] = deal(nan([mmax(ntr), matSz, nunits]));
sz = size(ct);

for i = 1:prod(matSz)
    % subscripts for 2D directions only
    % ii = diri==i & is2d;
    [dii, aii] = ind2sub(matSz, i);
    % Logical index of all stimulus presentations matching this direction
    ii = xyi==i;
    % Tediously ensure we don't mix up unit responses in this reshaping
    for u = 1:nunits
        ct(1:ntr(i), dii, aii, u) = rate.count(ii, u);
        tr(1:ntr(i), dii, aii, u) = rate.raw(ii, u);
    end
end

dv.tune.stimType = stimType;
dv.tune.rate = rate;
dv.tune.dirs = dirs;
dv.tune.spd = spd;
dv.tune.tr = tr;
dv.tune.trMu = squeeze(nanmedian(tr));
dv.tune.ctZ = squeeze(nanmean(ct))./squeeze(nanstd(ct));
dv.tune.commonParams = commonParams;
dv.tune.condPars = condPars;
dv.tune.condSet = condSet;


% % % ------------------------------------

%% Figure & Plotting params
fprintf('Plotting direction tuning')

msz = 8;
lw = 1.52;%1.8;

cols = get(0,'DefaultAxesColorOrder');
refcol = .91*[1 1 1];
xt = [-90:90:450];
xtl = compose('%0.5g', xt);
xtl(1:2:end) = {''};

% 'plotBoxAspectRatio'
pbaxr = [5,3,1];

% Subplot size & layout
spx = 2;
spy = nunits/spx;
% % !! %   spx == 2 makes subplots match stereo-probe config

% Double column (native stereo-probe config)
figW = 8;
figHscale = 3;
figsz = [figW, nunits/spx*figHscale];%nunits*figWscale, figHeight];

% subplot spacing ([a.u.])
spMargin = [.02 .02];


if ishandle(saveout)    % saveout>0 && ishandle(saveout)
    H1 = saveout;
    spbase = nunits;
elseif saveout<-1
    H1 = figureFS(abs(saveout), 'portrait', figsz);
    spbase = 0;
else
    H1 = figureFS([], 'portrait', figsz);
    spbase = 0;
end


if isempty( get(H1, 'name'))
    set(H1, 'name', [fileName,'_',stimType{2},'-bySpd'])
end


% Stimulus specifics
switch stimType{1}
    case 'gabors'
        % Only use 2d directions (...ignore twd & away conds for now)
        dd = dirs(:,1)==dirs(:,2);
        dirs2d = dirs(dd,1);
        
    case 'dotBall'
        dd = dirs(:,2)==0;
        dirs2d = dirs(dd,1);
        
end

% trMx = bsxfun(@rdivide, dv.rf.trMu, max(max(dv.rf.trMu )));
% rateCut = 2/3;

% dir range for tuning curve plots
if doPolar
    dirRng = 0:1:360;
    
else
    dirRng = -90:1:450; %-10:1:370;
    % wrapping index into dirs2d, dirs Xz, & data
    dlen = length(dirs2d);
    dirsII = [ find(dirs2d>270)', 1:dlen, find(dirs2d<90)']';
    dirsXax = [dirs2d(dirs2d>270)-360; dirs2d; dirs2d(dirs2d<90)+360];
    
    xl = [min(dirRng), max(dirRng)] + [-15,15]; %    xl = [min(dirRng), max(dirRng)+15];
end


%% Do the plotting
legendStr = {};
for u = 1:nunits
    figure(H1)
    titext = [];
    sp = subplot_tight(spy ,spx, u+spbase, spMargin, 'plotboxaspectratio',pbaxr);
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
    
        
    % xlabel(sprintf('ch.%d:\nO(%2.0f, %1.2f)\tD(%2.0f, %1.2f)', dv.uprb.id(u), [op,ob], [dp,db]), 'fontsize',10, 'verticalalignment','middle');
    xlabl = sprintf('ch.%d:\nO(%2.0f, %1.2f)    D(%2.0f, %1.2f)', dv.uprb.id(u), [op,ob], [dp,db]);
    %         titext = sprintf('ch.%d:\nO(%2.0f, %1.2f)\tD(%2.0f, %1.2f)\nO(%2.0f, %1.2f)\tD(%2.0f, %1.2f)', dv.uprb.id(u), [op,ob], [dp,db], [op3,ob3], [dp3,db3]);
    
    set(sp, 'linewidth',0.5)
    
    if doPolar
        set(sp, 'ThetaTick',thTic, 'ThetaTickLabel', thTicl, 'RTick', rticks, 'rlim',yl);
    else
        set(sp, 'plotBoxAspectRatio',pbaxr, 'XTick',xt,'XTickLabels',xtl,  'YTick', [rticks]);
        ylim(sp, yl); % [-15,375])
        xlim(sp, xl);   %[-15,375]);
    end

    
    dv.tune.ori(:,u) = [op, ob];
    dv.tune.dir(:,u) = [dp, db];
    
    if u==1
        try
            titext = sprintf('%s\t%s', fileName, dv.uprb.info.comment); %#ok<*AGROW>
        catch
            % old style
            titext = sprintf('%s\t%s', fileName, dv.uprb.info);
        end
        titext = {titext, sprintf('spd: %s deg/s', mat2str(spd))};
    end
    ht = title(titext, 'fontsize',10, 'fontweight','normal', 'interpreter','none');
    xlabel(xlabl, 'fontsize',10, 'fontweight','normal')
    
    % draw stim loc on rf plot
    try % prevent error crash if no associated rf plot
        stimPos = [commonParams.pos(1:2), commonParams.sz];
        ha = findobj(figure(rfFigH), 'type','axes','tag',sprintf('ch%drf',u));
        if ~isempty(ha)
            set(ha, 'NextPlot','add');  %axes(ha); hold on;
            hacl = get(ha,'clim');
            
            % Draw contours for gabor contrast envelope at [0.2, 0.5, 0.95] intensity
            gw = commonParams.sz; % stim fwhm    tstStim(3)/2.355;
            gww = gw*1.2; gx = linspace(-gww,gww,100); [gx,gy] = meshgrid(gx, gx); gg = zeros(size(gx)); gg(:) = hypot(gx(:), gy(:));   % gauss2D_R(gx, gy, gw, gw, 0, 1);
            contour(ha, gx+stimPos(1), gy+stimPos(2), gg, gw*[1 1], 'color',curveColor(end,:), 'linewidth',1);

            set(ha, 'clim',hacl);
            drawnow nocallbacks
        end
    end
    
    % Add waveform inset (to subplot[sp])
    % ...this should really be its own compartmentalized fxn.   2019: IT IS!!
    figure(H1)
    if ~isempty(dv.uprb.wf.mu)
        addWfInset(sp, dv.uprb.wf, u);
    end
    
    % progress
    fprintf('.')

end
fprintf('\n')


end %main function



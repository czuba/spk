function [dv] = tuneDispAwake(filebase, unitArgs, saveout, rfFigH)

% cd /Volumes/Tank1/projectData/adaptDist/awake/axel/data
figDir = '~/Dropbox/Science/projects/3dViewDist/kipp/figs';
    evalin('caller', sprintf('figDir = ''%s''', figDir))

%% Default inputs

% file name
if nargin<1 || isempty(filebase)
    filebase = []%
%     defarg('filebase', sprintf('%s_%s', 'kipp', datestr(now,'yyyymmdd')))
%     filebase = cell2mat(inputdlg('rfPos File Name String     (e.g. ''<subj>_YYYYMMDDa_<gridLoc>'')', 'FileName', 1, {filebase}));
end

% spikes to load/use
defarg('unitArgs', [666, 0]) % Default:  666=load all spikes, collapse all sorted & unsorted

if length(unitArgs)<2
    unitArgs = [unitArgs, 0];   % 2nd arg wf loading: 0=no[default], 1=yes & return ch stats,2=yes & return all wfs)
end

% save figs
defarg('saveout',nan)

if nargin<4 || isempty(rfFigH)
    rfFigH = [];
end

if ~exist('spkSrc','var')
    spkSrc = cellstr('uprb');
    %spkSrc = {'uprb'};
end

%% Load the files
if isstruct(filebase)
    % first input was existing dv struct, not filebase string
    dv = filebase;
else
    dv = initDaily(filebase, pwd, spkSrc);
    
    
    % Make sure this is data from the proper stimulus file type
    if isfield(dv.pds.baseParams.session, 'caller') && ~contains(lower(dv.pds.baseParams.session.caller.name), 'disptune')
        warning('PDS file generated from:  %s\nThis doesn''t appear to be the correct file type for computing Direction Tuning.', dv.pds.baseParams.session.caller.name)
        keyboard
    end
    
%     for i = 1:length(spkSrc)
%         dv.(spkSrc{i}) = plx_readerPar_Pldaps(dv.paths.(spkSrc{i}), unitArgs, dv.(spkSrc{i}).sync);
%         %dv.(spkSrc{i}) = getSpikes_plx(dv, spkSrc{i}, unitArgs);
%         
%         if isempty(dv.(spkSrc{i}).info)
%             dv.(spkSrc{i}).info = input(sprintf('Info string for this file (format:  [gridLoc]_[depth in microns] )\n\t:: '), 's');
%         end
%     end
end

fileName = dv.pds.baseParams.session.file(1:end-4);

% % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % 
% Following this point, only hacked up for "uprb" stim source
% % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % 


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
            commonParams.pos(2) = -commonParams.pos(2); % flip y-dimension when rendering texture stimuli
            %*** (...fix this in stim code to make universally 'correct': negative below horizontal meridian)
        commonParams.sz = dv.pds.baseParams.(baseModule).gabSz;
        commonParams.sf = dv.pds.baseParams.(baseModule).gabSf;
        %         condMat = cellfun(@(x) [(x.stimPos(1:2).*dv.pds.baseParams.(baseModule).gridSz(1:2) + dv.pds.baseParams.(baseModule).stimCtr(1:2)).*doYflip, x.dir] ...
        %             , dv.pds.condMatrix.conditions, 'uni',0);
        nDirs = rate.condDims(1);
        nSpd = rate.condDims(2);
        % Collapsing across speed [for quick & dirty]
        matSz = [nDirs, 1];

        % CONDITIONS:  [ LeftEyeDir, RightEyeDir, speedDegPerSec, posXpx, posYpx ]
        condMat = cellfun(@(x) [x.dir, x.gabTf ./ commonParams.sf, x.stimPos.*[1,-1] + commonParams.pos(1:2)], dv.pds.condMatrix.conditions, 'uni',0);
        condSet = cell2mat(condMat(:));
        
        % Stim conditions for ALL strobes (in order of presentation)
        strobeConds = condSet(rate.stimStrobes(:,2),:);
        
        % collapse across ori/direction
        if str2num(fileparts(pwd))>20181030
            [dirs, ~, diri] = unique(strobeConds(:,1:2), 'rows');
        else
            fprintf(2, '\n\nHACK FOR BAD BINO DIRECTION CONDITION!!!\n\n')
            [dirs, ~, diri] = unique(rate.condIdx(:,1), 'rows');
            dirs = condSet(1:nDirs, 1:2)
        end
        [xi, yi] = ind2sub(matSz, diri);
        
    case 'dotBall'
        commonParams.pos = dv.pds.baseParams.(baseModule).stimCtr;
        commonParams.gridSz = dv.pds.baseParams.(baseModule).gridSz;
        commonParams.sz = dv.pds.baseParams.(baseModule).ballSz;    %gabSz;
        commonParams.spd = dv.pds.baseParams.(baseModule).dotSpd;
        nDepth = rate.condDims(1);
        nDirs = rate.condDims(2);
        nSpd = length(commonParams.spd);   % 

        % Only frontoparallel directions at each disparity
        matSz = [nDepth, nDirs];

%         % CONDITIONS:  [ xyDir, motPlane, speedDegPerSec, posXpx, posYpx ]
%         condMat = cellfun(@(x) [x.dir, x.stimPos(1:2).*commonParams.gridSz(1:2) + commonParams.pos(1:2)], dv.pds.condMatrix.conditions, 'uni',0);
%         condSet = cell2mat(condMat(:));
        
        % get params shared across all modules, and ensure a 3rd (z) coord is present
        commonParams.pos = dv.pds.baseParams.(baseModule).stimCtr;
            if numel(commonParams.pos)<3, commonParams.pos(3) = 0; end
        commonParams.gridSz = dv.pds.baseParams.(baseModule).gridSz;
            %if numel(commonParams.gridSz)<3, 
                commonParams.gridSz(3) = 1; %end
        
        stimPos = cellfun(@(x) x.stimPos, dv.pds.condMatrix.conditions, 'uni',0);
        stimPos = cell2mat(stimPos(:)); % [nconds, xyz]
        if size(stimPos,2)<3, stimPos(:,3) = 0; end
        % combine condition position with grid and stim center
        stimPos = bsxfun(@plus, bsxfun(@times, stimPos, commonParams.gridSz), commonParams.pos);
        % NO y-axis flip needed for dots!!
        %         doYflip = [1,1, 1]; % flip y-dimension when rendering texture stimuli
        %         stimPos = bsxfun(@times, stimPos, doYflip);
        
        stimDir = cellfun(@(x) x.dir, dv.pds.condMatrix.conditions, 'uni',0);
        stimDir = cell2mat(stimDir(:)); % [nconds, motDirAxisRot]
        
        condSet = [stimPos, stimDir];
        
        
        % Stim conditions for ALL strobes (in order of presentation)
        strobeConds = condSet(rate.stimStrobes(:,2),:);
        
        % Convert strobe values into x & y indices, while collapsing across any
        % additional condMatrix dimensions (e.g. ori, direction, disparity, speed etc)

        % -- final [oi] output from ind2sub captures all [o]ther dimensions
        [dispi, moti, oi] = ind2sub(rate.condDims, rate.stimStrobes(:,2));
        [dmi] = sub2ind(matSz, dispi, moti); % 
%         is2d = axi==1;
        
        % Set of unique directions
        condPars = condSet(unique(dmi), 3:4);
        ds = unique(condPars(:,1), 'stable');
        ms = unique(condPars(:,2), 'stable');
%         ds = condPars(1:rate.condDims(1), 1);       % directions
%         ms = condPars(1:rate.condDims(1):end, 2);   % axes (??)
        

end


% trial count in each condition
ntr = reshape(hist(dmi, length(condPars)), matSz);
% ntr = reshape(hist(diri(is2d), length(dirs)), matSz);
% [ct tr] = deal(nan([mmax(ntr), nunits, matSz]));
[ct tr] = deal(nan([mmax(ntr), matSz, nunits]));
sz = size(ct);

for i = 1:length(condPars)
    % subscripts for 2D directions only
    % ii = diri==i & is2d;
    [dii, mii] = ind2sub(matSz, i);
    % Logical index of all stimulus presentations matching this direction
    ii = dmi==i;
    % Tediously ensure we don't mix up unit responses in this reshaping
    for u = 1:nunits
        ct(1:ntr(i), dii, mii, u) = rate.count(ii, u);
        tr(1:ntr(i), dii, mii, u) = rate.raw(ii, u);
    end
end

% % for i = 1:length(dirs)
% %     ii = diri==i;
% %     ct(1:ntr(i), :, i) = rate.count(ii,:);
% %     tr(1:ntr(i), :, i) = rate.raw(ii,:);
% % end
% % 
% % % move unit index to end for consistency/ease of use
% % ct = permute(ct, [1,3:length(sz),2]);
% % tr = permute(tr, [1,3:length(sz),2]);

dv.tune.stimType = stimType;
dv.tune.rate = rate;
dv.tune.condPars = condPars;
    dv.tune.depths = ds;
    dv.tune.disps = depth2disp(ds, dv.pds.baseParams.display.viewdist/100, dv.pds.baseParams.display.ipd)/60;
    dv.tune.dirs = ms;
dv.tune.tr = tr;
dv.tune.trMu = squeeze(nanmean(tr));
dv.tune.ctZ = squeeze(nanmean(ct))./squeeze(nanstd(ct));
dv.tune.commonParams = commonParams;

% % % ------------------------------------



%% Polar tuning plots
fprintf('Plotting Disparity tuning')

% trMx = bsxfun(@rdivide, dv.rf.trMu, max(max(dv.rf.trMu )));
% rateCut = 2/3;
figHeight = 8;
figWscale = 2;
figsz = [nunits*figWscale, figHeight];
spMargin = 1/nunits


if ishandle(saveout)    % saveout>0 && ishandle(saveout)
    H = figureFS(saveout, 'portrait', figsz);
    H2 = figureFS(saveout+1, 'portrait', figsz);
    H3 = figureFS(saveout+2, 'portrait', figsz);
    spbase = nunits;
elseif saveout<-1
    H = figureFS(abs(saveout), [], figsz);
    H2 = figureFS(abs(saveout)+1, [], figsz);
    H3 = figureFS(abs(saveout)+2, [], figsz);
    spbase = 0;    
else
    H = figureFS(400, [], figsz);
    H2 = figureFS(401, [], figsz);
    spbase = 0;
end
% subplot arrangements
spx = nunits;
if isempty(findobj(H, 'type','axes','tag',sprintf('ch%drf',1)))
    spbase = 0;
    spx = spx/2;
end
if isempty( get(H, 'name'))
    set(H, 'name', [fileName,'_mot'])
end
set(H2, 'name', [fileName,'_disp'])
set(H3, 'name', [fileName,'_dirDisp'])

switch stimType{1}
    case 'gabors'
        % Only use 2d directions (...ignore twd & away conds for now)
        dd = dirs(:,1)==dirs(:,2);
        dirs2d = dirs(dd,1);
        
    case 'dotBall'
        %         dd = dirs(:,2)==90;
        dirs2d = dv.tune.dirs;
        disps = dv.tune.disps;
        
        
end

legendStr = {};
for u = 1:nunits
% % % %     ii = ismember(xyi, find( trMx(:,:,u) >= rateCut ));
% % %     [db, dp, ob, op] = orivecfit(xy, dv.tune.trMu(:,u));
% % %     dv.tune.ori(:,u) = [op,ob];
% % %     dv.tune.dir(:,u) = [dp,db];
    
%     subplot_tight(4,3,u,.08);
%     subplot_tight(2,nunits, u+nunits, 0.3/nunits)
%     if nunits>1
%         % sp = subplot_tight(2,nunits, u+nunits, 0.3/nunits);
%         sp = subplot_tight(2,spx, u+spbase, 0.3/nunits);
%         %cla
%     else
%         sp = subplot_tight(2,spx, u+spbase, .1);
%         %cla
%     end
    
    % tuning limits
    yl = [-0.1, 1.1] *round2(mmax(dv.tune.trMu(:,:,u)), 5);
    if yl(2)<50, yl(2)=50; end  % set ylim >= 50 spk/s
    % vector sum tuning
    %     [db, dp, ob, op] = orivecfit(dirs(:,1), dv.tune.trMu(:,u));

    % Plot DIR TUNING for each disparity        (wrapping endpoint)
    xl = [-10,370];     xt = 0:90:360;
    % trTmp = dv.tune.trMu(dd,u);
    trTmp = dv.tune.trMu(:,:,u);
    figure(H);
    sp = subplot_tight(2,spx, u+spbase, spMargin);
    hold on;
    for D = 1:length(disps)
        [db(D), dp(D), ob(D), op(D)] = orivecfit(dirs2d, trTmp(D,:));
        ph = plot([0;dirs2d([1:end])], trTmp(D,[end,1:end]), '.-');
        %  plh = polarplot(d2r(dirsXz([end,1:end])), trTmp([end,1:end]),'.-');
    end
    hp = gca;
    ylim(yl)
    xlim(xl)
    curveColor = get(ph, 'color');
    box off
    set(hp, 'plotBoxAspectRatio', [2,3,1], 'xtick',xt);

    
    
    % Plot DISPARITY TUNING for each direction  (cartesian axes)
    figure(H2);
    xl = 1.1*disps([1,end]);        xt = -2:.5:2;
    sp2 = subplot_tight(2,spx, u+spbase, spMargin);
    hold on;
    for D = 1:length(dirs2d)
        %         set(sp, 'colororderindex', 2)
        ph = plot(disps, trTmp(:,D), '.-');
        %  plh = polarplot(d2r(dirsXz([end,1:end])), trTmp([end,1:end]),'.-');
    end
    hp = gca;
    ylim(yl);
    xlim(xl);
    curveColor = get(ph, 'color');
    box off
    set(hp, 'plotBoxAspectRatio', [2,3,1], 'xtick',xt);
        if u == 1,
            legendStr = {legendStr, 'disp'};
        end

        % %         % XZ tuning polar plot
        % %         plh = polarplot(d2r(dirs2d([end,1:end])), trTmp3([end,1:end]),'v-', 'markersize',2, 'linewidth',.5);
        % %         if u == 1,
        % %             legendStr = {legendStr, 'XZ'};
        % %         end
        % %
        % %         curveColor = get(plh, 'color');

        % rticks = unique(round2(max(trTmp)*[0.5,1], 5));
        % set(hp, 'ThetaTick',[0:45:360], 'RTick', rticks)
        
    % Plot DIR-DISP Surface         (wrapping dir endpoint)
    xl = [-10,370];     xt = 0:90:360;
    yl = 1.1*disps([1,end]);        yt = -2:.5:2;
    [xx, yy] = meshgrid([0;dirs2d(:)], disps);
    % trTmp = dv.tune.trMu(dd,u);
    trTmp = dv.tune.trMu(:,:,u);
    figure(H3);
    sp3 = subplot_tight(2, spx, u+spbase, spMargin);
    hold on;
    surf(xx, yy, trTmp(:,[end,1:end]), 'linestyle','none');%, 'FaceColor','interp');
    ylim(yl),   ylabel('disp', 'fontsize',8);
    xlim(xl),   xlabel('dir', 'fontsize',8);
    box off,    view(2)
    set(gca, 'plotBoxAspectRatio',[2,3,1], 'xtick',xt, 'ytick',yt);
    
    
    titext = sprintf('ch.%d:', dv.uprb.id(u));
%         titext = sprintf('ch.%d:\nO(%2.0f, %1.2f)\tD(%2.0f, %1.2f)', dv.uprb.id(u), [op,ob], [dp,db]);
%         titext = sprintf('ch.%d:\nO(%2.0f, %1.2f)\tD(%2.0f, %1.2f)\nO(%2.0f, %1.2f)\tD(%2.0f, %1.2f)', dv.uprb.id(u), [op,ob], [dp,db], [op3,ob3], [dp3,db3]);


% %             end
    set([sp,sp2,sp3], 'linewidth',0.5)
    
    dv.tune.ori(:,u) = [op, ob];
    dv.tune.dir(:,u) = [dp, db];
    
%     dv.tune.ori3(:,u) = [op3, ob3];
%     dv.tune.dir3(:,u) = [dp3, db3];
    
            %     ctS{u} =   sparse(yi, xi, rate.count(:,u),      matsz(1), matsz(2));
            %     durS{u} =  sparse(yi, xi, rate.trialdur,   matsz(1), matsz(2));
            %     trMuS{u} = sparse(yi, xi, rate.raw(:,u),        matsz(1), matsz(2));
            %     ntrS{u} =  sparse(yi, xi, 1,               matsz(1), matsz(2));
            
    if u==1
        titext = {sprintf('%s\t%s ', fileName, dv.uprb.info.comment), titext}; %#ok<*AGROW>
        % ht = title(sprintf('%s\t%s', fileName, dv.uprb.info), 'fontsize',10, 'interpreter','none');
%         if nunits>=4
%             set(ht,'horizontalalignment','left')
%         else, set(ht,'horizontalalignment','center'), end
    end
    ht = title(sp, titext, 'fontsize',10, 'fontweight','normal', 'interpreter','none');
    ht = title(sp2, titext, 'fontsize',10, 'fontweight','normal', 'interpreter','none');
    ht = title(sp3, titext, 'fontsize',10, 'fontweight','normal', 'interpreter','none');
    
    
    %% draw stim loc on rf plot
    try % prevent error crash if no associated rf plot
        stimPos = [commonParams.pos(1:2), commonParams.sz];
        ha = findobj(figure(rfFigH), 'type','axes','tag',sprintf('ch%drf',u));
        set(ha, 'NextPlot','add');  %axes(ha); hold on;
        hacl = get(ha,'clim');
        % tune stimulus (conv. dia to radius)
%         circle(stimPos(1), stimPos(2), stimPos(3)/2, 'none', .8*[1 1 1], 2)
        tstCol = .8*[1 .4 .4];
        %           circle(tstStim(1), tstStim(2), tstStim(3)/2.355/2, 'none', tstCol, 1)
        %            tstStim = [7.5, -0.5, 3.5*2.355]
%         % Draw contours for gabor contrast envelope at [0.2, 0.5, 0.95] intensity
%         gw = commonParams.sz; % stim fwhm    tstStim(3)/2.355;
%         gww = gw*2.355; gx = linspace(-gww,gww,100); [gx,gy] = meshgrid(gx, gx); gg = gauss2D_R(gx, gy, gw, gw, 0, 1);
%         contour(ha, gx+stimPos(1), gy+stimPos(2), gg, [.2 .5 .95], 'color', tstCol);

        % Draw contours for gabor contrast envelope at [0.2, 0.5, 0.95] intensity
        gw = commonParams.sz; % stim fwhm    tstStim(3)/2.355;
        gww = gw*1.2; gx = linspace(-gww,gww,100); [gx,gy] = meshgrid(gx, gx); gg = zeros(size(gx)); gg(:) = hypot(gx(:), gy(:));   % gauss2D_R(gx, gy, gw, gw, 0, 1);
        contour(ha, gx+stimPos(1), gy+stimPos(2), gg, gw/2*[1 1], 'color',curveColor, 'linewidth',1);
        % set(gco, 'linewidth',1);
        
        %         uistack(ha, 'bottom')
        set(ha, 'clim',hacl);
        drawnow nocallbacks
    end
    
    figure(H2);
    % add waveform inset
    if ~isempty(dv.uprb.wf.ci)
        ap = get(hp, 'position');
        aw = axes('position', [ap(1)+0.8*ap(3), ap(2), 0.4*ap(3), 0.2*ap(4)]);
        plot(aw, dv.uprb.wf.ci(:,:,u), 'r', 'linewidth',.5);
        hold on,
        plot(aw, dv.uprb.wf.mu(:,u), 'k', 'linewidth',.5);        
        if range(dv.uprb.wf.mu(:))<1
            % wf vals in mV
            yl = .1; %yl = .04;
        else
            % else...wf vals in arb. a/d units (3000 is good)
            yl = 3000;
        end   
        yt = yl/2;
        set(aw, 'color','none', 'yAxisLocation','right', 'linewidth',0.5, 'xticklabel',[],...
            'yticklabel',[], 'ytick',yt*[-1,0,1], 'ylim',yl*[-1,1],...
            'xlim',[1,size(dv.uprb.wf.ci,1)], 'box','off', 'plotboxaspectRatio',[3,2,1]);
        uistack(aw, 'top')
    end
    drawnow nocallbacks
    % progress
    fprintf('.')
end
fprintf('\n')
% saveFigTriplet(1, get(gcf,'name'), 0)

% dv.rf.d2pix = GetConversionFactor(dv.expo, 'deg', 'pix');%full(dv.expo.environment.Conversion.units.ConversionMatrix(6,4));
% xdata = round(expodata(:,5)*dv.rf.d2pix);
% ydata = round(expodata(:,6)*dv.rf.d2pix);


%% run hacky plotting script
% rfPlotHak(dv, h, saveout)



return
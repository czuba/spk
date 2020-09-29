function [dv] = tuneXzAwake(filebase, unitArgs, saveout, rfFigH)


if evalin('caller', sprintf('exist(''figDir'',''var'') && ~isempty(figDir)'))
    figDir = evalin('caller', 'figDir')     %#ok<*NOPRT>
end
% % Push [figDir] upstream
% figDir = '~/Dropbox/Science/projects/3dViewDist/kipp/figs';
%     evalin('caller', sprintf('figDir = ''%s''', figDir))

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

if saveout>0 || ishandle(saveout)
    rfFigH = saveout;
elseif nargin<4 || isempty(rfFigH)
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
    dv = initDaily(filebase, pwd, spkSrc, unitArgs);
    
    
    % Make sure this is data from the proper stimulus file type
    if isfield(dv.pds.baseParams.session, 'caller') && ~contains(lower(dv.pds.baseParams.session.caller.name), 'xztune')
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
ntr = reshape(hist(xyi, length(dirs)), matSz);
% ntr = reshape(hist(diri(is2d), length(dirs)), matSz);
% [ct tr] = deal(nan([mmax(ntr), nunits, matSz]));
[ct tr] = deal(nan([mmax(ntr), matSz, nunits]));
sz = size(ct);

for i = 1:length(dirs)
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
dv.tune.dirs = dirs;
dv.tune.tr = tr;
dv.tune.trMu = squeeze(nanmean(tr));
dv.tune.ctZ = squeeze(nanmean(ct))./squeeze(nanstd(ct));
dv.tune.commonParams = commonParams;

% % % ------------------------------------



%% Polar tuning plots
fprintf('Plotting polar tuning')

if ishandle(saveout)    % saveout>0 && ishandle(saveout)
    H = saveout;
    figure(H); % make active
    spbase = nunits;
else
    H = figure;
    spbase = 0;
end

if isempty( get(H, 'name'))
    set(H, 'name', [fileName,'_tune'])
end

% subplot arrangements
spx = nunits;
if isempty(findobj(H, 'type','axes','tag',sprintf('ch%drf',1)))
    spbase = 0;
    spx = spx/2;
end

if spbase==0 && nunits~=32
    spx = 8;
    spy = ceil(nunits/spx);
else
    spx = 16;
    spy = 2;
end

% plotbox aspect ratio
pbax = [1 1 1]; %for polar plots  [3,2,1];

% % figHeight = 8;
% % figWscale = 2;
% % figsz = [nunits*figWscale, figHeight];

figHeight = 4*spy;
figWidth = 2+3*spx;
figsz = [figWidth, figHeight];

spMargin = [.03 .025];   %1/nunits

% apply size formatting to [each] figure
figureFS(H, 'portrait', figsz);


switch stimType{1}
    case 'gabors'
        % Only use 2d directions (...ignore twd & away conds for now)
        dd = dirs(:,1)==dirs(:,2);
        dirs2d = dirs(dd,1);
        
    case 'dotBall'
        dd = dirs(:,2)==90;
        dirs3d = dirs(dirs(:,2)==90, 1);
        dirs2d = dirs(dirs(:,2)==0,1);
        respIdx = 2;
end

lcols = get(groot, 'DefaultAxesColorOrder');


legendStr = {};
for u = 1:nunits

    figure(H);

    if nunits>1
        % sp = subplot_tight(2,nunits, u+nunits, 0.3/nunits);
        sp = subplot_tight(spy, spx, u+spbase, spMargin);
        %cla
    else
        sp = subplot_tight(spy, spx, u+spbase, .1);
        %cla
    end
    
    % vector sum tuning
%     [db, dp, ob, op] = orivecfit(dirs(:,1), dv.tune.trMu(:,u));
    % trTmp = dv.tune.trMu(dd,u);
    if ndims(dv.tune.trMu)==3
        trTmp2 = dv.tune.trMu(:, 1, u);
        trTmp3 = dv.tune.trMu(:, 2, u);
    else
        %trTmp = dv.tune.trMu(dirs(:,2)==90, u);
        trTmp2 = dv.tune.trMu(dirs(:,2)==0, u);
        trTmp3 = dv.tune.trMu(dirs(:,2)==90, u);
    end
    
    
    [db2, dp2, ob2, op2] = orivecfit(dirs3d, trTmp2);
    [db3, dp3, ob3, op3] = orivecfit(dirs3d, trTmp3);

    yl = [-0.1, 1.1] *round2(mmax([trTmp2(:);trTmp3(:)]), 5);

        % Plot 2D (fronto) tuning in BLUE on polar axes, wrapping endpoint        % set(gca, 'colororderindex', 2)
        if ~isempty(dirs2d)
            plh = polarplot(d2r(dirs2d([end,1:end])), trTmp2([end,1:end]),'.-');
            hp = gca; hold on;
            set(plh, 'color', lcols(3,:), 'linewidth',1, 'markersize',8);
            
            curveColor = get(plh, 'color');
            
            if u == 1,
                legendStr{end+1} = '2d';
            end
        end
        
        if ~isempty(dirs3d)
            % Plot XZ tuning in ORANGE on polar axes, wrapping endpoint        % set(gca, 'colororderindex', 3)
            plh = polarplot(d2r(dirs3d([end,1:end])), trTmp3([end,1:end]),'.-');
            hp = gca; hold on;
            set(plh, 'color', lcols(4,:), 'linewidth',1, 'markersize',8);
            rticks = unique(round2(max(trTmp3)*[0.5,1], 5));
            set(hp, 'ThetaTick',[0:45:360], 'RTick', rticks)
            
            curveColor = get(plh, 'color');
            
            if u == 1,
                legendStr{end+1} = '3d';
            end
        end
        rticks = unique(round2(max([trTmp2(:); trTmp3(:)])*[0.5,1], 5));
        set(hp, 'ThetaTick',[0:45:360], 'RTick', rticks)        

        if u == nunits
            lh = legend(legendStr, 'Location','south','orientation','horizontal','FontSize',10,'LineWidth',2);
            lhpos = get(lh,'position');
            set(lh, 'position',[.85,.01,lhpos(3:4)]);
        end
        
%         ph = plot([0;dirsXz([1:end])], trTmp([end,1:end]), '.-');
%         hp = gca; hold on;
%         ylim(yl)
%         curveColor = get(ph, 'color');
%         box off
%         set(hp, 'plotBoxAspectRatio', [2,3,1], 'xtick',[0:45:360]);
%         xlim([-15,375])


        % %         % XZ tuning polar plot
        % %         plh = polarplot(d2r(dirs2d([end,1:end])), trTmp3([end,1:end]),'v-', 'markersize',2, 'linewidth',.5);
        % %         if u == 1,
        % %             legendStr = {legendStr, 'XZ'};
        % %         end
        % %
        % %         curveColor = get(plh, 'color');

        
        % xlabel(sprintf('ch.%d:\nO(%2.0f, %1.2f)\tD(%2.0f, %1.2f)', dv.uprb.id(u), [op,ob], [dp,db]), 'fontsize',10, 'verticalalignment','middle');
        %         titext = sprintf('ch.%d:\nO(%2.0f,%1.2f) D(%2.0f,%1.2f)', dv.uprb.id(u), [op,ob], [dp,db]);
        if mod(dv.uprb.id(u),1) == 0
            % old unit id format: chan*10+unit#
            titext = sprintf('ch.%d',dv.uprb.id(u));
        else
            % kilosort id:  Unit#.peakChannel#
            titext = sprintf('u%05.2f',dv.uprb.id(u));
        end
        % add tuning details
        titext = [titext, sprintf(':\n2d(%2.0f,%1.2f)  ||  3d(%2.0f,%1.2f)', [dp2,db2], [dp3,db3])];
%         titext = sprintf('ch.%d:\nO(%2.0f, %1.2f)\tD(%2.0f, %1.2f)\nO(%2.0f, %1.2f)\tD(%2.0f, %1.2f)', dv.uprb.id(u), [op,ob], [dp,db], [op3,ob3], [dp3,db3]);


% %             end
    set(hp, 'linewidth',0.5)
    
    dv.tune.ori2d(:,u) = [op2, ob2];
    dv.tune.dir2d(:,u) = [dp2, db2];
    dv.tune.ori3d(:,u) = [op3, ob3];
    dv.tune.dir3d(:,u) = [dp3, db3];
    

    if u==1
        titext = {fileName, [dv.tune.stimType{2},' || ',dv.uprb.info.comment], titext}; %#ok<*AGROW>
    end
    
    ht = title(titext, 'fontsize',10, 'fontweight','normal', 'interpreter','none');
        if nunits>=2 && u==1
            set(ht,'horizontalalignment','left')
        else,
            set(ht,'horizontalalignment','center'),
        end

        
    % draw stim loc on rf plot
    try % prevent error crash if no associated rf plot
        stimPos = [commonParams.pos(1:2), commonParams.sz];
        ha = findobj(figure(rfFigH), 'type','axes','tag',sprintf('ch%drf',u));
        set(ha, 'NextPlot','add');  %axes(ha); hold on;
        hacl = get(ha,'clim');
        tstCol = .8*[1 .4 .4];

        % Draw contours for gabor contrast envelope at [0.2, 0.5, 0.95] intensity
        gw = commonParams.sz; % stim fwhm    tstStim(3)/2.355;
        gww = gw*1.2; gx = linspace(-gww,gww,100); [gx,gy] = meshgrid(gx, gx); gg = zeros(size(gx)); gg(:) = hypot(gx(:), gy(:));   % gauss2D_R(gx, gy, gw, gw, 0, 1);
        contour(ha, gx+stimPos(1), gy+stimPos(2), gg, gw/2*[1 1], 'color',curveColor, 'linewidth',1.5, 'linestyle','--');
        
        %         uistack(ha, 'bottom')
        set(ha, 'clim',hacl);
        drawnow nocallbacks
    end
    
    figure(H);
    % add waveform inset
    if ~isempty(dv.uprb.wf.mu)
        addWfInset(hp, dv.uprb.wf, u);
    end
    
    drawnow nocallbacks
    % progress
    fprintf('.')
end
fprintf('\n')
% saveFigTriplet(1, get(gcf,'name'), 0)


return

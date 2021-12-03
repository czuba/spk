function [dv, Hout] = rfPosAwake3d(baseName, unitArgs, addto, doFitting, varargin)
% function [dv, Hout] = rfPosAwake3d(baseName, unitArgs, addto, doFitting, ['param1',val1, 'param2',val2, ...])
% 
% Inputs are getting super gangly...
% 
% INPUTS:
%   [baseName]  =   String containing a [unique] portion of PDS file name you'd like to process
%                   -- the last 4-digit number in PLDAPS file name (experiment start time) is good
%                   -- Also accepts a cell of {name, path} for calling from another location (i.e. pwd ~= recording session directory)
%   [unitArgs]  =   Flags controlling how/what spike data is loaded (see below)
%   [addto]     =   Figure number to use: if exist, will add to that figure, if not, will create new figure with that number
%   [doFitting] =   Fit circular gaussians [true]
%   PV-pairs    =   Parameter value pairs [...not yet implemented]
% 
% OUTPUTS:
%   [dv] struct is the primary data variables output struct for all Czuba analysis code.
%   -- This struct syncs and pairs together data paths, stimulus, spike, and info structures in a standardized Czuba-format
%   -- initDaily.m is the catchall function for initializing a dv struct
% 
% 
% Suggested usage:
%   dvRf = rfPosAwake3d('1527', [-666, 1], -1)
%   -- '1527'       is the timestamp of PDS file interested in (e.g. <file>_1527.PDS)
%   -- [-666, 1]    [-666] will collapse all units, negative forces use of unsorted spike source file
%                   [1]    will load waveform means & ci for each unit
%   -- -1           plotting flag: negative triggers new figure to be used; positive will use that figure number, if it already
%                   exists it will attempt to add new plots to existing axes based on axis tags (expert mode; YMMV)
%                   addto>0 mostly for plotting multiple rfPosition files in same figure (e.g. different eye or region)
% 
% 
%   [unitArgs] flag options  (...hackily complex, but is what it is)
%       default: [666, 0]  Loads all spikes (including unsorted), no waveforms;
%       arg(1) == 0:    All sorted & unsorted spikes
%       arg(1) == 1:    Ignore unsorted spikes (don't load PLX sortcode == 0)
%       arg(1) == 666:  Collapse all spike sort codes  (default; good for ignoring junk units created during recording)
%       arg(2) == 0:    Dont load waveforms. (default)
%       arg(2) == 1:    Load waveforms, but only pass out wf mean & ci of each unit (default)
%       arg(2) == 2:    Load waveforms, pass out wf mean, ci, & every individual wf (this will be big!)
%       arg(2) == 3;    Load waveforms, pass out wf mean & means from 4 temporal intervals across file
%                       (good for coarse stability check)
% 
% See also: initDaily, syncPlexon2PDS, calcRate_condMatrix, figureFS,
%           plx_readerPar_Pldaps, getKiloPath, getSpikes_kilo
% 
% 201x-xx-xx  TBC  Evolved over years from acute to awake experiments, and beyond
% 2020-01-29  TBC  Cleaning & commenting. Add record of calling string in output figure.UserData (see end of code)
% 

%  TODO:
%   parameter-value pairs & defaults:
%   'spx'       = n subplot columns     [2]
%   'spy'       = n subplot rows        [nunits/spx]
%   


%% Default inputs
% fit RFs to individual or averaged trial responses?
useTrialAverages = 1;
% plot interpolant surface, or basic response image?
plotInterpolant = 1; 


if evalin('caller', sprintf('exist(''figDir'',''var'') && ~isempty(figDir)'))
    figDir = evalin('caller', 'figDir')     %#ok<*NOPRT>
end

% file name
if nargin<1 || isempty(baseName)
    baseName = [] 
    basePath = pwd;

elseif iscell(baseName)
    % more specific inputs to initDaily
    switch length(baseName)
        case 1
            basePath = pwd;
        case 2
            basePath = baseName{2};
    end
    baseName = baseName{1};
else
    basePath = pwd;
end

%** This relies on defarg.m; an ancient piece of code for checking/setting default values.
%   Its "just don't look" code that isn't broke, so... --TBC 2020-01

% spikes to load/use
defarg('unitArgs', [666, 0]); % Default:  666=load all spikes, collapse all sorted & unsorted

% add rfs to already existing fig?
defarg('addto', []);

if nargin<4 || isempty(doFitting)
    doFitting = 1;
elseif doFitting && isstruct(baseName)
    % remove previous tunning & fit outputs
    baseName = rmfield(baseName,'rf');
elseif ~doFitting && ~isstruct(baseName)
    % Incompatible inputs...tuning/fits needed
    fprintf(2, '~!~\tFitting not requested, but inputs do not appear sufficient to skip it.\n\tReverting to default doFitting=1\n')
    doFitting = 1;
end

% TODO:  Parse additional PV-pair inputs

if ~exist('spkSrc','var')
    spkSrc = cellstr('uprb');
    %spkSrc = {'uprb'};
end

%% Load the files:  create [dv] struct
% % [dv] struct is the primary data variables output struct for all Czuba analysis code.
% % -- This struct syncs and pairs together data paths, stimulus, spike, and info structures in a standardized Czuba-format
% % -- initDaily.m is the catchall function for initializing a dv struct
% % 
% % Typically a session will consist of multiple [dv] structs, one for each RF/tuning analysis component.
% % -- In .mat data files, these are delineated by suffix on the 'dv', as in:  dvRf, dvT[une], dvD[isparity],...etc
% % -- dv_output.mat files are a single data file containing multiple dv structs,
% %    and ideally the output table from >> tt = scanSesh; containing info about each PLDAPS file in that session
% % 
if isstruct(baseName)
    % first input was existing dv struct, not baseName string
    dv = baseName;
    clear baseName % free up memory from dupes
else
    dv = initDaily(baseName, basePath, spkSrc, unitArgs);
end

if isfield(dv.pds.baseParams.session, 'caller') && ~contains(lower(dv.pds.baseParams.session.caller.name), 'rfpos')
    warning('PDS file generated from:  %s\nThis doesn''t appear to be the correct file type for computing RF position.', dv.pds.baseParams.session.caller.name)
    keyboard
end

fileName = dv.pds.baseParams.session.file(1:end-4);


% % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % 
% Following this point, only coded for PLDAPS condMatrix & "uprb" stim source
% % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % 

%% Compute rate from condMatrix

% calc spike rate for each presentation
rate = calcRate_condMatrix(dv.uprb, dv.pds);

nunits = size(rate.count,2);

%% Parse condMatrix conditions & collapse across nessary dimensions
% % Dependent on structure of ACTUAL MATRIX MODULE    (...in this case:  glDraw.doRfPos_gabGrid.m)
% baseModule = dv.pds.condMatrix.modNames.currentStim{1};

% Type of RF mapping stimulus used (e.g. 'gabors' or 'dotBall')
stimType = dv.info.stimType;
baseModule = stimType{1};

% ...here we go.

xyzSize = rate.condDims(1:3);

switch stimType{1}
    case 'gabors'
        doYflip = [1,-1]; % flip y-dimension when rendering texture stimuli
        %*** (...fix this in stim code to make universally 'correct': negative below horizontal meridian)
        % ...done, I think, but hardly ever use gabor rf pos anymore, since 3Ddots so much better
        condMat = cellfun(@(x) [(x.stimPos(1:2).*dv.pds.baseParams.(baseModule).gridSz(1:2) + dv.pds.baseParams.(baseModule).stimCtr(1:2)).*doYflip, x.dir] ...
            , dv.pds.condMatrix.conditions, 'uni',0);
        
    case 'dotBall'
        % flag for depth already in visual degrees
        depthDeg = isfield(dv.pds.condMatrix.conditions{1}, 'stimDepthDeg');

        stimCtr = dv.pds.baseParams.(baseModule).stimCtr;
        gridSz = dv.pds.baseParams.(baseModule).gridSz;
        viewDist = dv.pds.baseParams.display.viewdist;
        ipd = dv.pds.baseParams.display.ipd;
        
        if depthDeg % isfield(dv.pds.condMatrix.conditions{1}, 'stimDepthDeg')
            % updated version with disparity offset defined in deg horizontal disparity
            condMat = cellfun(@(x) [([x.stimPos(1:2).*gridSz(1:2), viewDist + disp2depth(60*x.stimDepthDeg, viewDist/100, ipd)] + stimCtr), x.dir], dv.pds.condMatrix.conditions, 'uni',0);
        else
            % %         if ndims(dv.pds.condMatrix.conditions)>3 % has [x,y,z,dir,...]
            % NOTE: Z is not transformed by gridSz, instead must pass through and add as offset (stimCtr too) to viewing distance.
            condMat = cellfun(@(x) [([x.stimPos(1:2).*gridSz(1:2), x.stimPos(3)+viewDist]  + stimCtr), x.dir] ...
                , dv.pds.condMatrix.conditions, 'uni',0);
        end
        % %         else
        % %             condMat = cellfun(@(x) [(x.stimPos .* dgridSz + stimCtr) .* doYflip, x.dir] ...
        % %                 , dv.pds.condMatrix.conditions, 'uni',0);
        % %         end
        %         condMat = cellfun(@(x) [(x.stimPos(1:2).*dv.pds.baseParams.(baseModule).gridSz(1:2) + dv.pds.baseParams.(baseModule).stimCtr(1:2)).*doYflip, x.dir] ...
        %             , dv.pds.condMatrix.conditions, 'uni',0);
    otherwise
        fprintf(2, '~~~\tError!!! Unrecognized RF stimulus type:\t%s\n\tCrash imminent...\n', stimType{1});
        keyboard
end
condSet = cell2mat(condMat(:));

% Stim conditions for ALL strobes (in order of presentation)
strobeConds = condSet(rate.stimStrobes(:,2),:);


% Convert strobe values into x & y indices, while collapsing across any
% additional condMatrix dimensions (e.g. ori, direction, disparity, etc)

% -- final [~] output from ind2sub captures all additional dimensions
[xi, yi, zi, oi] = ind2sub(rate.condDims, rate.stimStrobes(:,2));
[xyzi] = sub2ind([xyzSize], xi, yi, zi);
% Set of unique XY positions properly sorted
xyz = condSet(unique(xyzi), 1:3);
xs = unique(xyz(:,1),'stable'); %xyz(1:rate.condDims(1),1);
ys = unique(xyz(:,2),'stable'); %xyz(1:rate.condDims(1):end,2);
zs = unique(xyz(:,3),'stable');
% NOTE: unique(strobeConds(:,1:2), 'rows')  returns an improperly sorted unique set(!)

% trial count in each condition
ntr = reshape(hist(xyzi, length(xyz)), xyzSize);
[ct tr] = deal(nan([mmax(ntr), xyzSize, nunits]));
sz = size(ct);

for i = 1:length(xyz)
    % subscripts for this xyz index
    [xii, yii, zii] = ind2sub(xyzSize, i);
    % Logical index of all stimulus presentations matching this xy position
    ii = xyzi==i;
    % Tediously ensure we don't mix up unit responses in this reshaping
    for u = 1:nunits
        ct(1:ntr(i), xii, yii, zii, u) = rate.count(ii,u);
        tr(1:ntr(i), xii, yii, zii, u) = rate.raw(ii,u);
    end
end

dv.rf.stimType = stimType;
dv.rf.rate = rate;
dv.rf.xyz = xyz;
    dv.rf.xs = xs;
    dv.rf.ys = ys;
    dv.rf.zs = zs;
dv.rf.tr = tr;

dv.rf.trMu = squeeze(nanmedian(tr));
dv.rf.trVar = squeeze(nanvar(tr));
dv.rf.ctZ = squeeze(nanmean(ct))./squeeze(nanstd(ct));

% % % % 
%% TODO:
% Ideally should have a confirmation step 1.1 to ensure averaging only happens
% within fully consistent conditions.
%       (e.g. condMatrix location defined relative to stimCtr, but stimCtr was following
%       the mouse position!)
% 
% confirm all xy pos & directions *actually* presented from .data{}
% --Sync is already linked to trial in 3rd column of sync.strobes,
% but how do we tie syncs to module identities? (e.g. 'gabors01', 'gabors02',...)
% --Cant just reshape, b/c syncs need not be from completed trials, so partial
% trials would break things.
% --stimStrobes already has other syncs removed, just add 4th column with counting
% index for each repetition of trial # in 3rd column.
%   --modName = sprintf('%s%02d', dv.pds.condMatrix.modNames.currentStim, stimStrobes(:,4));
% 
% % % % 


%% Polar plots of responsivity
% Make or find the appropriate figure

% Subplot size & layout
[spx, spy, figsz, spMargin, pbaxr] = mkProbeSubplots(nunits);
%   mkProbeSubplots == scrappy attempt to standardize probe data subplotting
%   - optional second input flag:
%       'g' (default) produces 8-by-n grid of subplots
%       'v' for 2-by-n vertical pairs; arranged like unsorted channels of stereoprobe
%       'h' for n-by-2 rows; marginally more 'monitor-friendly', primarily used back in the <=16 channel days)

figOri = 'landscape';

% Crufty flag trying to do too much...[new figure, choose number, existing figure,...lots of edge cases]
if isempty(addto) || ~ishandle(addto)
    if ~isempty(addto) && isnumeric(addto)
        % Use [addto] as figure number
        if addto<0 && ishandle(abs(addto))
            close(abs(addto));
        end
        H1 = figureFS(abs(addto), figOri, figsz);
    else
        % Open new figure
        H1 = figureFS([], figOri, figsz);
    end
    % name it
    set(H1, 'name', [fileName,'_rf'], 'tag',figDir)
else
    % 3rd argin is handle to existing figure
    H1 = addto;
end


%%  Fit 2d rfs
fprintf('\nFitting RFs...')

trRaw = num2cell(rate.raw, 1);

% normalize trial means to peak response (add eps to prevent divide by zero)
trMx = bsxfun(@rdivide, dv.rf.trMu, max(max(dv.rf.trMu ))+eps);
trDims = ndims(dv.rf.trMu);

% Handles grid sizes with or without depth planes!!
iDims = 1:trDims-1;
% trial stim positions
trPos = strobeConds(:, iDims);  % resp.pos

% mean tr matrix  % maybe not used?...
imU = num2cell(dv.rf.trMu, iDims);
trMu = squeeze(num2cell(dv.rf.trMu, iDims));
xyzMat = dv.rf.xyz;
trVar = squeeze(num2cell(dv.rf.trVar, iDims));
trVarMu = cellfun(@(x) nanmean(x(:)), trVar);

%   mmxy = trMx(:,:,u);
mmxy = squeeze(num2cell(trMx, iDims))';

if trDims<4
    % dummy flags
    ud = dv.info.viewDist;
    uii = ones(length(trPos));
else
    % unique depths
    [ud, ~, uii] = unique(trPos(:,3));
    mmxy = cellfun(@(x) sum(x,3), mmxy, 'uni',0);
end


fxn = @gaussian3d_1; % 3d circular gaussian (forced sigma-x == sigma-y)
rfOpts = optimset('LargeScale','on','display','off');

if useTrialAverages
    trXY = num2cell(xyzMat(:,1:2),1); % XY pos for trial averages
else
    trXY = num2cell(trPos(:,1:2), 1); % XY pos for individual trials
end


% initialize fit outputs
[rfFit] = deal(cell(nunits,1));

tic
parfor u = 1:nunits,
                    
    % index of max xy
    [x0i, y0i] = find(mmxy{u}==max(mmxy{u}(:)));    
    if length(x0i)>1
        x0i = x0i(fix(length(x0i)/2));
        y0i = y0i(fix(length(y0i)/2));
    end
    % max xy pos
    x0 = xs(x0i);
    y0 = ys(y0i);
    
    % init sigma is eccentricity
    s0 = hypot(x0,y0);
    spc = min([diff(xs);diff(ys)]);
    
    if useTrialAverages
        yvals = trMu{u}(:); % use trial medians
    else
        yvals = trRaw{u}; % use individual trial rates
    end
    
    % initial estimate:  [baseline, amp, x-mean, y-mean, sigma]     ...if 2d var allowed: [...x-sigma, y-sigma]
    fit0 = [mean(yvals), max(yvals)-mean(yvals),  x0,   y0, s0];
    % bounds
    LB = [0,    0,  xs(1)-spc,  ys(1)-spc,  s0/4];
    UB = [mean(yvals)*2, max(yvals)*2,    xs(end)+spc,    ys(end)+spc, s0*4];
    
    % Fit XY Gaussian RF!
    [out resnorm] = lsqcurvefit(fxn, fit0, trXY, yvals, LB, UB, rfOpts);
    
    % fit quality
    r2 = 1-resnorm / ((length(yvals)-1)*var(yvals));
    % Also equiv to:  1-resnorm/sum((trRaw{u}-mean(trRaw{u})).^2)
    % ...but not good for nonlinear fits, esp not with multiple depth params.

    % package outputs
    rfFit{u} = [out, r2, resnorm];
    
    % progress
    fprintf('.')
end
sprintf('done. (%3.3f sec/unit)\n\n', toc/nunits)

% Unpack parfor cell outputs
dv.rf.fxn = fxn;
dv.rf.fxnPars = {'base', 'amp', 'xmu', 'ymu', 'sigma'};
dv.rf.fit = cell2mat(cellfun(@(x) x(1:end-2), rfFit, 'uni',0));
dv.rf.r2 = cellfun(@(x) x(end-1), rfFit);
dv.rf.resnorm = cellfun(@(x) x(end), rfFit);


%% Fit quality
% Computing this metric is jacked because:
%   --fit is non-linear
%   --isn't optimized for any single depth condition
%   --variance for non-responsive depths/conds is particularly low (...reliably zero!)
% 
%   % % Print out rf fits in command window:  [unit, x, y, ecct, fwhm, r2]
%   rftable = rfPosFitTable(dvRf.rf)


%% Plot RFs
fprintf('\nPlotting RFs... ')

% mapping grid lims
gl = [min(xyz); max(xyz)];
% make image pixels fall in corrrect locations
gl = gl + kron(diff(gl)./xyzSize./2, [-1;1]); 

% plotting axes limits
rfStimLims = gl(1:4) + 3*[-1,1,-1,1]; % make stimulus edges visible
rfAxLims = [-5,30,-30,5];  % [-8,32,-25,10]; % a standard large area of screen

    % Set axes lims inclusive of stim & standard
    ii = [1,3];     rfAxLims(ii) = min([rfStimLims(ii); rfAxLims(ii)]);
    ii = [2,4];     rfAxLims(ii) = max([rfStimLims(ii); rfAxLims(ii)]);
    
    % manual override, zoom on stimulus area
%     rfAxLims = rfStimLims;

% xy grid for plotting fits
spc = min([diff(xs);diff(ys)]);
[xxRF, yyRF] = meshgrid(xs(1)-10:.2:xs(end)+10, ys(1)-10:.2:ys(end)+10);


figure(H1)

% find any pre-existing axis tags
ha = findobj(H1, 'type','axes');

for u = 1:nunits,

    if isempty(addto) || addto<0 || isempty(ha)
        sp = subplot_tight(spy, spx, u, spMargin);
        % tag as channel RF
        set(sp, 'tag',sprintf('ch%drf',u), 'NextPlot','add');
    else
        % Retrieve tagged channel rf axes & hold on
        sp = ha(strcmp(get(ha, 'tag'), sprintf('ch%drf',u)), 'NextPlot','add');
    end
            
    % Response matrix for this unit
    if ndims(dv.rf.trMu)==3
        % just xy
        imU = dv.rf.trMu(:,:, u)';
    elseif ndims(dv.rf.trMu)==4
        % xyz, collapse z
        imU = dv.rf.trMu(:,:, :, u);
    else
        error('rfPosAwake3:badDimensions','Unrecognized RF response matrix dimensions:\t%s\n\t%s', mat2str(size(dv.rf.trMu)), dv.paths.pds)
    end
    
    
    %% Plot it
    rez = 3; % upsampling

    if plotInterpolant
        % Use matlab Interpolant fit for generating RF image/surface
        warnbak = warning('off', 'MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId');
        rez = ceil(rez*length(dv.rf.xs));
        
        if ndims(dv.rf.trMu)==3     % RFs sampled in single plane (presumably in plane of fixation)
            % Fit data directly
            V = dv.rf.trMu(:,:,u);    %dv.rf.trMu(:,:,:,u);
            rawXY = strobeConds(:,1:2);
            F = scatteredInterpolant(rawXY(:,1), rawXY(:,2), dv.rf.rate.raw(:,u), 'natural');
            % Get interpolated response surface
            [xx, yy] = ndgrid(linspace(gl(1,1), gl(2,1),rez), linspace(gl(1,2), gl(2,2),rez)); % zz = 0*ones(size(xx));
            Vi = F(xx, yy);
            hh = surf(xx, yy, -Vi, 'linestyle','none', 'FaceColor','interp');
            view(2); axis tight equal, axis(rfAxLims)

            % Flip colormap for negative spike rate plot (so overlays show above surf plot)
            colormap(gca, flipud(parula(64))); 
            set(gca, 'clim', -[max(Vi(:)), min(Vi(:))])
            
            %             % make flat view of RF surface from below so that overlaid reference elements are visible
            %             view(0,-90); axis tight equal, axis(rfAxLims), set(gca, 'ydir','reverse')
            
        elseif ndims(dv.rf.trMu)==4     % RFs sampled in multiple depth planes
            % Fit data directly
            V = dv.rf.trMu(:,:,:,u);
            rawPos = strobeConds(:,1:3);
            [ud, ~, uii] = unique(rawPos(:,3));
            
            % If multiple depth planes, plot each in separate figure/subplot
            if length(ud)>1
                if ~exist('rfTabs', 'var')
                    rfTabs = mkFigureTabs( (1:nunits)+800, dv.info.pdsName);
                end
                % common color range
                icl = prctile(dv.rf.rate.raw(:,u), [5,95], 'all');
                if diff(icl)<1, icl = [0,1]; end
                
                figureFS(rfTabs(u), 'portrait',[6,4]);  % (u+100);
                set(rfTabs(u), 'name',sprintf('%s-%02d',get(H1,'name'),u));
                di = [];
                for d = 1:length(ud)
                    % just this depth
                    ii = uii==d;
                    %  F = scatteredInterpolant(rawPos(ii,1), rawPos(ii,2), dv.rf.rate.raw(ii,u), 'natural');
                    rrd = trMu{u}(:,:,d); rrdn = numel(rrd);
                    F = griddedInterpolant( reshape(xyzMat(1:rrdn,1),size(rrd)), reshape(xyzMat(1:rrdn,2),size(rrd)), rrd, 'makima');
                    % Get interpolated response surface
                    [xx, yy] = ndgrid(linspace(gl(1,1), gl(2,1), rez), linspace(gl(1,2), gl(2,2),rez)); % zz = 0*ones(size(xx));
                    Vi = F(xx, yy);
                    
                    di(d) = subplot(1,length(ud), d);
                    hh = surf(xx, yy, Vi, 'linestyle','none', 'FaceColor','interp');
                    view(2); axis tight equal, % axis(rfAxLims)
                    title(sprintf('u%d VD=%2.2fcm',u, dv.rf.zs(d)));
                end
                %             pause(.01) % pause is faster than processing lag (???)
                %             linkprop(di,{'CLim','View'});  %...doesn't work inside loop(?!?)
                icl = cell2mat(get(di, 'clim'));
                icl = [min(icl(:,1)), max(icl(:,2))];
                set(di, 'clim',icl);
                colorbar(di(2),'southoutside')
            end
            % back to the main fig window
            figure(H1);
            % ignore disparity and interp all as though same depth
            F = scatteredInterpolant(rawPos(:,1), rawPos(:,2), dv.rf.rate.raw(:,u), 'natural');

            % Get interpolated response surface
            [xx, yy] = ndgrid(linspace(gl(1,1), gl(2,1),rez), linspace(gl(1,2), gl(2,2),rez)); % zz = 0*ones(size(xx));
            Vi = F(xx, yy);
            hold on
            hh = surf(xx, yy, -Vi, 'linestyle','none', 'FaceColor','interp');
            view(2); axis tight equal, axis(rfAxLims)
                        
            % Flip colormap for negative spike rate plot (so overlays show above surf plot)
            colormap(gca, flipud(parula(64))); 
            set(gca, 'clim', -[max(Vi(:)), min(Vi(:))])
            
            %             % flip Y and view from bottom so xyunity ref lines are visible (...maybe more trouble than its worth)
            %             view(0,-90); axis tight equal, axis(rfAxLims), set(gca, 'ydir','reverse')

        end
        warning(warnbak);

    else
        % Plot rf image map
        if ndims(imU)==2
            imU = imresize(imU, rez, 'method','lanczos2');
            imsz = size(imU);
            
            imagesc(gl(:,1), gl(:,2), imU, 'parent',sp);
            set(sp, 'ydir','normal', 'xtick',[-60:10:60], 'ytick',[-40:10:40], 'fontsize',10);
            axis equal;  axis(sp, rfAxLims);
        else
            
            imU = imresize(imU, rez, 'method','lanczos2');
            imsz = size(imU);
            if length(dv.rf.zs)>1
                % plot each depth plane in separate figure/subplot  (...slow & hacky)
                if ~exist('rfTabs', 'var')
                    rfTabs = mkFigureTabs( (1:nunits)+100, dv.info.pdsName);
                end
                figure(rfTabs(u));
                for i = 1:length(dv.rf.zs)
                    di(i) = subplot_tight(1, length(dv.rf.zs), i, .07);
                    imagesc(gl(:,1), gl(:,2), imU(:,:,i)');
                    title(sprintf('u%d VD=%2.2fcm',u,dv.rf.zs(i)));
                end
                set(di, 'ydir','normal', 'fontsize',10);
                axis equal tight; box off;
                linkaxes(di); linkprop(di, {'clim'});
            end
            imU = max(imU,[],3)';
            imagesc(gl(:,1), gl(:,2), imU, 'parent',sp);     %imagesc(xy(:,1), xy(:,2), imU);
            figure(H1)
            set(sp, 'ydir','normal', 'xtick',[-60:5:60], 'ytick',[-40:5:40], 'fontsize',10);
            axis equal;  axis(sp, rfAxLims);
        end
    end
    xyunity; box off

    % Plot rfFit contour
    if isfield(dv.rf, 'fxn')
        fout = dv.rf.fxn(dv.rf.fit(u,:), {xxRF(:),yyRF(:)});
        rf50 = dv.rf.fit(u,1) + dv.rf.fit(u,2)/2;
        contour(xxRF, yyRF, reshape(fout, size(xxRF)), [1,1]*rf50, 'color','k','linestyle',':', 'xliminclude','off', 'yliminclude','off');
    end

    xlabel(sp, sprintf('ch %.2f\n(%4.1f, %4.1f), %4.1f',  dv.uprb.id(u), dv.rf.fit(u,3:4), 2.355*sqrt(dv.rf.fit(u,5))), 'fontsize',10);

    if u==1
        % full title on first plot only
        titext = {sprintf('%s  VD:%d ', dv.info.pdsName, dv.info.viewDist), [dv.info.stimType{2},' || ',dv.uprb.info.comment]}; %#ok<*AGROW>};
        ht = title(sp, titext, 'fontsize',10, 'interpreter','none');
        if nunits>=2
            tialign = 'left';
        else
            tialign = 'center';
        end
        set(ht,'horizontalalignment','center')

    elseif u==nunits
        % add stim location grid to last image
        plot(sp, dv.rf.xyz(:,1), dv.rf.xyz(:,2), '.','color',.97*[1 1 1], 'markersize',2)
    end
    
    
    % Add waveform inset (to subplot[sp])
    % ...this should really be its own compartmentalized fxn.   2019: IT IS!!
    figure(H1)
    if ~isempty(dv.uprb.wf.mu)
        addWfInset(sp, dv.uprb.wf, u, [0,-.1]); % shift inset up slightly
    end
    
    % progress
    fprintf('.')
end
fprintf('\n')
    
drawnow

% ensure focus on RF plot
axes(sp);

Hagg = scratch_plotRf_agg(dv);

if nargout>1
    % output figure handles
    Hout = [H1, Hagg];
    if exist('rfTabs', 'var')
        Hout = [Hout, rfTabs];
    end
end


% Record calling command in UserData
% HACKY:  only command line calls for now
dbSt = dbstack;
if numel(dbSt)==1
    try
        % Get session history from Java object, convert to cell string:
        historyText = cellstr(char(com.mathworks.mlservices.MLCommandHistoryServices.getSessionHistory));
        % Extract last string:
        callingStr = historyText{end};
        % Only continue if contains this file name
        if contains(callingStr, mfilename)
            set(H1, 'UserData', callingStr)
        end
    end
end



% % % % % % % % % %
%% Nested Functions for gangly code snippets
% % % % % % % % % %




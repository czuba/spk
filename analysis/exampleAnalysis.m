function [dv, Hout] = exampleAnalysis(baseName, unitArgs, addto, doFitting, varargin)
% function [dv, Hout] = exampleAnalysis(baseName, unitArgs, addto, doFitting, ['param1',val1, 'param2',val2, ...])
% 
% Grepped from rfPosAwake3d.m; inputs are getting super gangly, but initialization portion of code will get
% spike & stimulus data loaded & synced, and show how to retrieve spike times for all/specific conditions
% 
% INPUTS:
%   [baseName]  =   String containing a [unique] portion of PDS file name you'd like to process;
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
% For example usage:
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
%       arg(2) == 1:    Load waveforms, but only pass out wf mean & ci of each unit
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


%% Parse inputs & set defaults

% Check if figure directory is already defined in calling workspace
% -- relevant to saveFigTriplet.m
if evalin('caller', sprintf('exist(''figDir'',''var'') && ~isempty(figDir)'))
    figDir = evalin('caller', 'figDir')     %#ok<*NOPRT>
else
    figDir = [];
end

% Best to return handles of figures for saving, rather than write EVERY analysis code to handle saving separately
% See saveFigTriplet.m for a good figure saving utility
Hout = [];


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

% Update this check with the expected PLDAPS calling function name
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

% Using knowledge of the conditions matrix from the stimulus code & dv.pds.condMatrix
% (see each dv.pds.condMatrix.conditions for condition-specific parameter settings),
% parse the rate struct outputs into conditions, average across trials,
% refine windows, etc
% 
% NOTE:  Use particular caution in reshaping condMatrix indices into conditions. Dimension
%   transpositions/confusions can be subtle but damning in your analysis. 


%% Match up condition numbers with stimulus parameters of interest
% Compute final parameters AS PRESENTED to the subject by compiling shared & unique stimulus parameters of each condition
% 
% % !!! % %
% NOTE: unique(strobeConds(:,1:2), 'rows')  returns an improperly sorted unique set(!)
% ... MUST use unique( ..., 'stable') on each row
% % !!! % % % % %

% -- RF pos XYZ example:
% % %         xyzSize = rate.condDims(1:3);
% % % 
% % %         % Get shared stimulus parameters are from dv.pds.baseParams:
% % %         stimCtr = dv.pds.baseParams.(baseModule).stimCtr;
% % %         gridSz = dv.pds.baseParams.(baseModule).gridSz;
% % %         viewDist = dv.pds.baseParams.display.viewdist;
% % % 
% % %         % *** Tricky/convoluted cellfun call ***
% % %         %  - XY [stimPos] of each condition are stimulus grid offsets from a shared center postion [stimCtr]
% % %         %  - But Z is not transformed by gridSz, instead must pass through and add as offset (stimCtr too) to viewing distance.
% % %         condMat = cellfun(@(x) [([x.stimPos(1:2).*gridSz(1:2), x.stimPos(3)+viewDist]  + stimCtr), x.dir] ...
% % %                 , dv.pds.condMatrix.conditions, 'uni',0);
% % % 
% % %         condSet = cell2mat(condMat(:));
% % % 
% % %         % Stim conditions for ALL strobes (in order of presentation)
% % %         strobeConds = condSet(rate.stimStrobes(:,2),:);
% % % 
% % % 
% % %         % Convert strobe values into x & y indices, while collapsing across any
% % %         % additional condMatrix dimensions (e.g. ori, direction, disparity, etc)
% % % 
% % %         % -- final [~] output from ind2sub captures all additional dimensions
% % %         [xi, yi, zi, oi] = ind2sub(rate.condDims, rate.stimStrobes(:,2));
% % %         [xyzi] = sub2ind([xyzSize], xi, yi, zi);
% % % 
% % %         % Set of unique XY positions, properly sorted
% % %         xyz = condSet(unique(xyzi), 1:3);
% % %         xs = unique(xyz(:,1),'stable'); %xyz(1:rate.condDims(1),1);
% % %         ys = unique(xyz(:,2),'stable'); %xyz(1:rate.condDims(1):end,2);
% % %         zs = unique(xyz(:,3),'stable');
% % % 
% % %         % trial count in each condition
% % %         ntr = reshape(hist(xyzi, length(xyz)), xyzSize);
% % %         [ct tr] = deal(nan([mmax(ntr), xyzSize, nunits]));
% % %         sz = size(ct);
% % % 
% % %         for i = 1:length(xyz)
% % %             % subscripts for this xyz index
% % %             [xii, yii, zii] = ind2sub(xyzSize, i);
% % %             % Logical index of all stimulus presentations matching this xy position
% % %             ii = xyzi==i;
% % %             % Tediously ensure we don't mix up unit responses in this reshaping
% % %             for u = 1:nunits
% % %                 ct(1:ntr(i), xii, yii, zii, u) = rate.count(ii,u);
% % %                 tr(1:ntr(i), xii, yii, zii, u) = rate.raw(ii,u);
% % %             end
% % %         end


% Compile the outputs
dv.examp.stimType.stimType = stimType;
dv.examp.rate = rate;


% - Here is the further processing of above RF example:
% dv.rf.stimType = stimType;
% dv.rf.rate = rate;
% dv.rf.xyz = xyz;
%     dv.rf.xs = xs;
%     dv.rf.ys = ys;
%     dv.rf.zs = zs;
% dv.rf.tr = tr;
% 
% dv.rf.trMu = squeeze(median(tr,'omitnan'));
% dv.rf.trVar = squeeze(var(tr,'omitnan'));
% dv.rf.ctZ = squeeze(mean(ct,'omitnan'))./squeeze(std(ct,'omitnan'));




% % % % 
%% Stimulus specific Analysis code
% 
% Continue analyses on compiled spike counts, rates, and means/medians
% See rfPosAwake3d.m or tuneDispAwake.m for examples of further analysis
% I don't currently have a good example of PSTH or STA style timecourse
% analysis code using this pipeline, but hopefully the approaches taken
% in the examples above are good fodder for development...
%   -- TBC 2020-01
% 



% % % % % % % % % % % %
% % % % % % % % % % % %
%% Tuning fit convention
% % % % % % % % % % % %
% % % % % % % % % % % %
% Very useful convention of placing fitted function handle in output struct

% For example:
fxn = @gaussian3d_1; % 3d circular gaussian (forced sigma-x == sigma-y)
% fitOpts = optimset('LargeScale','on','display','off');
% [out resnorm] = lsqcurvefit(fxn, fit0, xdat, ydat, LB, UB, fitOpts);
dummyfit = [30, 150, 2.5, -13, 25;   10, 80, 6, -8, 12];

% Fitted outputs with function
dv.examp.fxn = fxn;
dv.examp.fxnPars = {'base', 'amp', 'xmu', 'ymu', 'sigma'};
dv.examp.fit = dummyfit; % nUnit-by-nParams fitted outputs

if 0
    % Easily plot fitted results using the function handle
    u = 1;
    [xx,yy] = meshgrid(-5:.5:20, 5:-.5:-20);
    
    fitout = dv.examp.fxn( dv.examp.fit(u,:), {xx,yy});
    figure;
    surf(xx,yy,fitout);
    view(2); axis image;
end
% % % % % % % % % % % %
% % % % % % % % % % % %




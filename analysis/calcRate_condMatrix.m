function [rate] = calcRate_condMatrix(spkSrc, pds, late, syncID, cm)  %, baseWin, expo)
% function [rate] = calcRate_condMatrix(spkSrc, pds, late, syncID)  %, baseWin, expo)
% 
% Compute spike [rate] struct on a given [spkSrc] struct (i.e. dv.uprb) based on
% the [pds] PLDAPS stimulus struct.
% 
% Optional inputs:
%   [late]      latency windows (can be series of n-by-2 start & end times relative to sync)
%               defaults to stimulus duration
%   [syncID]    Best unused for now...use fields of rate struct ot parse after-the-fact
% 
% OUTPUTS:
%   [rate]  standardized czuba-struct of 
% 
% 
% NOTE:  Syncs aren't sent till a matrix stimulus OFFSET!
%   This is by design (...for better or worse) because that way only completed presentations produce a sync.
%   However, this means response windows must be definied relative to the OFFSET TIMESTAMP.
% 
% 
% 2020-01-29  TBC  Commenting. Needs clean-up, but 'just works' for now
% 


if ~exist('late','var') || isempty(late)
    % DEFAULT: full duration that matrix stimulus module(s) were active
    %   shifted by respLatency of spkSrc (def = 50 ms)
    late = []; % latency (in sec)
end

if ~exist('syncID','var') || isempty(syncID)
    % DEFAULT: computes rate on all condMatrix presentations
    % Optional logical index for subset of syncs used to compute rate
    % ...best not to mess with this. Just parse/combine conditions after rate struct is returned
    syncID = [];
%     % % % elseif islogical(syncID)
%     % % %     syncID = find(syncID);
end

% allow separate condMatrix input (for multi-condMatrix experiments)
if nargin<5 || isempty(cm)
    cm = pds.condMatrix;
else
    if numel(cm)>1
        warning('Must input a single condMatrix for rate calc')
        keyboard
    end
end

if isfield(cm,'indexRange')
    % upper limit of strobed condition indices: cm.baseIndex+(0:cm.indexRange)
    indexRange = cm.indexRange;
else
    indexRange = numel(cm.conditions);
end

%% Get spike data

syncs = spkSrc.sync.strobe(:,1);
rate.id  = spkSrc.id;
rate.snr    = spkSrc.snr;


%% Identify syncs from condMatrix
baseIndex = pds.condMatrix.baseIndex;
stimPass = spkSrc.sync.strobe(:,2)>baseIndex & spkSrc.sync.strobe(:,2)<=(baseIndex+indexRange);
if ~isempty(syncID)
    if islogical(syncID)
        % syncID is logical subset of all strobes
        stimPass = stimPass & syncID;
    else
        % syncID is subset of strobe **values** (i.e. NOT index numbers)
        stimPass = stimPass & ismember(spkSrc.sync.strobe(:,2)-baseIndex, syncID);
    end
end
    
% unused, but passed out in [rate] struct
stimStrobes = spkSrc.sync.strobe(stimPass, :);
stimStrobes(:,2) = stimStrobes(:,2)-baseIndex;

    
% condMatrix stim duration IN SECONDS
tdur = double(diff(pds.baseParams.(pds.condMatrix.modNames.matrixModule{1}).modOnDur));
    % This looks convoluted, but just takes diff of onset & offset time from first matrixModule.
    % --!!-- Assumes duration of all modules are the same...

% Response latency IN SECONDS
if isfield(spkSrc,'respLatency')
   respLatency = spkSrc.respLatency; 
else
    % dummy estimate for latency of response
    % ...hardcoded, but at least recorded in rate struct:-/
    respLatency = 50/1000;
end


% % % % % 
% WEIRDNESS AHEAD:  Syncs aren't sent till a matrix stimulus OFFSET!
% This is by design (...for better or worse) because that way only completed presentations produce a sync.
% However, this means response windows must be definied relative to the OFFSET TIMESTAMP.
% % % % % 
if isempty(late)
    % use full stimulus duration by default
    late = [-tdur, 0];
end

late = late + respLatency;

% % calc spike rate for each presentation
% rate = calcRate(spkSrc, late, stimPass);

%% compute spike rate (parallelized)
fprintf('Computing spike rate in [%3.0f:%3.0f]ms window...\n', late'*1000);

% slice vars for parallel processing
tNum = rate.id;

spikeT = spkSrc.ts;     % spikeT is a cell [nunits, 1] of spike timestamps for each unit/chan.

% ensure class agreement (rare hitch due to mixed single/double timing types)
late = cast(late, 'like', spikeT{1});

% Parse requested syncs..
if ~isempty(stimPass)
    if islogical(stimPass)
        if length(stimPass)~=length(syncs)
            error('Whoa...length of logical stimPass not equal to number of syncs!\nThis should be a logical index.\tCrashing...')
        end
        sID = find(stimPass);
    elseif max(stimPass)<=length(syncs)
        sID = stimPass;
    else
        error('stimPass inputs to %s are incorrect!\nRequested indices greater than number of spike syncs present.', mfilename)
    end
else
    % compute rate on ALL condMatrix stimulus presentations     [DEFAULT]
    sID = find(stimPass);
end

nlate = size(late,1);
[rC, rT] = deal(cell(length(spikeT), nlate));
[rC{:}, rT{:}] = deal(nan(length(sID),1));

tic
fprintf('\n')

%% Parallel loop
parfor t = 1:length(spikeT)
    fprintf('%3.0f, ',t);%disp(t)
    for si = 1:length(sID)
        for l = 1:nlate
            s = sID(si);
            rC{t,l}(si,1) = sum( spikeT{t} >= (syncs(s)+late(l,1))  &  spikeT{t} <= (syncs(s)+late(l,2)) );              %#ok<PFBNS>
            rT{t,l}(si,1) = rC{t,l}(si,1)  /  (diff(late(l,:)));
        end
    end
    
end

fprintf('\n\n%3.2f sec per unit\n\n',toc/size(tNum,1))
ns = length(sID);
sz = size(rT);
% reshape parfor cells to matrices
rate.count  = reshape( cell2mat(rC), [ns, sz]);
rate.raw    = reshape( cell2mat(rT), [ns, sz]);


rate.trialdur = diff(late'); %trialdur;
rate.respLatency = respLatency;
rate.late = late;

%% New stuff for condMatrix
% nunits = size(rate.count,2);

% Convert condMatrix indices of each sync into multi-dim subscripts
% ...here we go.
rate.stimStrobes = stimStrobes;
rate.nStrobes = size(stimStrobes,1);
rate.condDims = size(pds.condMatrix.conditions); % size of original matrix needed to interpret strobeVals (i.e. don't resize this to match any syncID subselection)

% strobe indexes to condMatrix subscripts
condIdx = cell(size(rate.condDims));
[condIdx{:}] = ind2sub(rate.condDims, stimStrobes(:,2));
% convert to matrix (...cell needed to capture all outputs of ind2sub)
rate.condIdx = cell2mat(condIdx);


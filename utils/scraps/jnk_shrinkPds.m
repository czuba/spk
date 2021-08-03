
% load all pds fields
tic
pds = load(origName, '-mat');%, 'baseParams','data');
toc

structfun(@(x) getByteSize(x, 'mb'), pds)

%%

if isfield(pds,'conditions')
    nc = length(pds.conditions);
    if isfield(pds.conditions{nc-1},'display')
        % remove large [.display] field from [virtually unused] pds.conditions struct
        % - is all in either pds.baseParams.display or pds.data{}.display already
        pds.conditions = cellfun(@(x) rmfield(x,'display'), pds.conditions, 'uni',0);
    end
end

%%
% prune excess fields from 'dotBall' stimuli (RFpos, 3D tuning; ...pds files >100MB)
fn = pds.baseParams.pldaps.modNames.all(contains(pds.baseParams.pldaps.modNames.all, {'ddDots', 'xzDots'}));
if ~isempty(fn)
    for i = 1:length(pds.data)
        for ii = fn
            if isfield(pds.data{i}.(ii{:}),'pos0')
                pds.data{i}.(ii{:}) = rmfield(pds.data{i}.(ii{:}),'pos0');
            end
            if isfield(pds.data{i}.(ii{:}),'pos')
                pds.data{i}.(ii{:}) = rmfield(pds.data{i}.(ii{:}),'pos');
            end
        end
    end
end


%%
disp('pds byte size')
disp('-------------')
structfun(@(x) getByteSize(x, 'mb'), pds)
disp('-------------')
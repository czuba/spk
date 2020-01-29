function kilopath = getKiloPath(basePath, varargin)
% function kilopath = getKiloPath(basePath, varargin)
% 
% Find matching kilosort path, given basePath to a raw spike file.
%   Example:  basePath = dv.paths.(spkSrc{ss})
% 
% [optional second input with alternate target directory name]
% 
% 2019-02-21  TBC  Wrote it.


[spkPath, spkFile, spkExt] = fileparts(basePath);

% allow substitute (i.e. merged sort)
if nargin>1
    spkFile = varargin{1};
    if isstring(spkFile)
        spkFile = char(spkFile);
    end
end

% Search for match
kilopath = [];  tries = 0;
while isempty(kilopath) 
    ks = fullfile(spkPath, 'KiloSort', ['*',spkFile(1:end-tries),'*']);
    kilopath = dir( ks );
    tries = tries+1;
    if (length(spkFile)-tries)<3
        warning('No KiloSort directories found for:  %s\n\t...crash likely', ks)
        kilopath = [];
        return
    end
end

% There can be only one...
if length(kilopath)>1
    [kilopath, ok] = chooseFile(kilopath, 'Multiple matching KiloSort directories found. Select one:', 'single');
    if ~ok
        warning('A %d kilosort directories found, but user canceled or selected none.\n\t...returning empty', length(kilopath))
        kilopath = [];
        return
    end
end

% Return full path to most recent match
[~,ii] = max([kilopath.datenum]);

kilopath = fullfile(kilopath(ii).folder, kilopath(ii).name);

end %main function
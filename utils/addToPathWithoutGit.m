function err = addToPathWithoutGit(dir, excludes, withSubdirs)
% function addToPathWithoutGit(dir, excludes, withSubdirs)
%
% Add [dir] to Matlab path without including all the hidden versioning junk (.git, .svn)
% -- Optionally also exclude subdirs matching any entries in [excludes] string (or cell of strings)
% -- by default includes all subdirectories, set [withSubdirs] to false to only add the main [dir]
%
% 2021-07-14  TBC  Modernization & Added optional error status output
% 

if nargin<3 || withSubdirs
    a = genpath(dir);
    withSubdirs = ' and subdirectories';
else
    a = dir;
    withSubdirs = [];
end

if isempty(a)
    fprintf('%s not found...attempting to continue\n', dir);
    if nargout>0
        err = true;
    end
    return
else
    b=textscan(a,'%s','delimiter',':');
    b=b{1};
    b(contains(b, '.git')) = [];
    b(contains(b, '.svn')) = [];
    if nargin>1
        if ~iscell(excludes), excludes = {excludes}; end
        for i = 1:numel(excludes)
            if ~isempty(excludes{i})
                b(contains(b, excludes{i})) = [];
            end
        end
    end
    addpath(b{:})
    disp([dir, withSubdirs, ' added to the path']);
    if nargout>0
        err = false;
    end
end

end % main function


function [isSym, nfo] = isSymlink(fn)
% helper function to determine if input file/directory name is a symlink
%
% 2021-12-15  tbc  Wrote it.
% 


if ischar(fn)
    fn = cellstr(fn);
elseif isstruct(fn) && isfield(fn,'name')
        % standard output of dir()
        fn = {fn.name};
elseif ~iscell(fn)
    keyboard % unknown input format
end
n = numel(fn);

isSym = false(1,n);

if IsWindows
    fprintf(2, '\n!!!\t%s has not been coded for Windows use \n!!!\tBy default, will return "false" & attempt to continue...\n\n', mfilename)

elseif isunix
    for i = 1:n
        [err, statMsg] = system(sprintf('stat %s', fn{i}));
        if ~err
            isSym(i) = contains(statMsg, 'symbolic link', 'IgnoreCase',true);
        else
            fprintf(2, statMsg)
            isSym(i) = false;
        end

        if nargout>1
            nfo{i} = statMsg; %#ok<AGROW>
        end
    end
end

end %main function

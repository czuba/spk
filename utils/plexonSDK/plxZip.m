function splits = plxZip(plxdir, base, outname, exclude)
% function splits = plxZip(plxdir, base, outname, exclude)
% 
% Concatenates all plx files with <base> in their file name (or all .plx files in plxdir if
% base is not provided) and saves as 'outname' (outname defaults to the last dir in plxdir).
% 
% Returns split points in 'splits' (for use with plx_split.c), and also saves the splits and orig
% filenames as outname_splits.mat for use with plxUnzip.m (the counterpart of this fxn)
% 
% 10-16-2012  TBC  Wrote it.

% Check inputs & set defaults
if ~exist('plxdir','var') || isempty(base)
    plxdir = pwd;
end

if ~exist('base','var') || isempty(base)
    base = [];
end

if ~exist('outname','var') || isempty(outname)
    if isempty(base)
        [~,outname] = fileparts(plxdir);
        outname = [outname,'_zip'];
    else
        outname = [base,'_zip.plx'];
    end
end

if ~strcmpi(outname(end-3:end),'.plx')
    outname = [outname,'.plx'];
end


%% tell the world your story
fprintf([repmat('-',1,80),'\nProcessing plx files in %s\n'], plxdir)

%% get list of source plx files in plxdir
srcf = dir([plxdir,filesep,'*',base,'*.plx']);

% Exclude certain plx files based on string in expo file name
%       ...reading this block of code may lower your IQ
if exist('exclude','var')
    if ~iscell(exclude), exclude = {exclude}; end % cellify
    fprintf('\nWARNING: plx files must have paired expo xml file to be detected for exclusion!!\n')
    fprintf('\tExcluding ''%s'' file(s):\t',exclude{:})
    exponames = dir([plxdir,'/*.xml']);
    clear badegg kill
    for e = 1:size(exponames)
            badegg(e) = ~isempty(regexp(exponames(e).name, sprintf('%s|',exclude{:}), 'once'));
            % badegg(e) = ~isempty(strfind(exponames(e).name, exclude));
    end
    badegg = {exponames(badegg).name};
    % check each file for matching run value
    kill = false(size(srcf));
    for s = 1:length(srcf)
        for e = 1:length(badegg)
            if strfind(srcf(s).name,['p',badegg{e}(strfind(badegg{e},'#')+1:strfind(badegg{e},'[')-1)])
                kill(s) = true;
                fprintf('%s,\t',srcf(s).name);
            end
        end
    end
    fprintf('\n\n')
    % trim any hits
    srcf(kill) = [];
end

% List files
fprintf('\t%s\n', srcf.name)

%% Smash them together
% ...a little dirty b/c we have to make and delete partials as we go...
splits = zeros(1,length(srcf)-1);

if length(srcf)<3
    % there's only one pair to combine...
    fprintf('\tMerging...\n\t\t%s & %s\n',srcf(1).name, srcf(2).name);
    splits(1) = plx_merge(srcf(1).name, srcf(2).name, outname);
else
    % the first pair
    fprintf('\tMerging...\n\t\t%s & %s\n',srcf(1).name, srcf(2).name);
    thatout = ['part1-',outname];
    splits(1) = plx_merge(srcf(1).name, srcf(2).name, thatout);
    
    for i = 3:length(srcf)
        if i~=length(srcf)
            % all the ones inbetween
            fprintf('\tMerging...\n\t\t%s & %s\n', srcf(i).name, thatout)
            thisout = ['part',num2str(i-1),'-',outname];
            splits(i-1) = plx_merge(thatout, srcf(i).name, thisout);
            delete(thatout);
            thatout = thisout;  % Seuss Code....
        else
            % the last one and clean up
            fprintf('\tMerging...\n\t\t%s & %s...last one\n', srcf(i).name, thatout)
            splits(i-1) = plx_merge(thatout, srcf(i).name, outname);
            delete(thatout);
        end
    end
end

%% Save what you did so you can undo it later
zipout = [outname(1:end-4),'_splits.mat'];
save(fullfile(plxdir,zipout), 'splits', 'srcf', 'outname');
fprintf(['\nCombined plx saved at:  %s\n\tZip:\t%s\n\tSplits:\t%s\n',repmat('-',1,80),'\n'], plxdir, outname, zipout)
        
    
function plxUnzip(plxdir, splitfile, zipname, savedir)
% function plxUnzip(plxdir, splitfile, zipname, savedir)
% Deconstructs a merged plx file(s) based on the split mat file created by plxZip.m
% (outname_splits.mat...see plxZip help)
% 
% Places resulting files in ./sort directory (creates if doesn't exist) and appends
% current date to new file names so that you don't accidentally overwrite the originals
% 
% 10-16-2012  TBC  Wrote it.

% Check inputs & set defaults

try
    load(fullfile(plxdir,splitfile))
catch
    [splitfile, plxdir] = uigetfile('*splits.mat', 'Select *_split.mat file (...saved out by plxZip)');
    if isempty(splitfile)
        error('Must provide splitfile (output from plxZip.m) to direct unzip action...')
    end
    load(fullfile(plxdir,splitfile))
end

if ~exist('zipname','var') || isempty(zipname)
    [zipname, plxdir] = uigetfile([outname(1:end-4),'*.plx'], 'Select [sorted] .plx file to split into original parts:');
end

% output directory
if ~exist('savedir','var')
    savedir = './sort/';
end
if ~exist(savedir,'dir')
    mkdir(savedir)
end

% sorted component file string (...wish was shorter, but this is unique and informative)
fend = ['_sort',datestr(now,'yymmdd'),'.plx'];  % just incase this fxn spans midnight...

% tell the world your story
fprintf('plxZip file:\t%s\n\tSplit into...\n', fullfile(plxdir,zipname))

% slice up the zip into the original components
% ...a little dirty b/c we have to make and delete partials as we go...
if length(splits)==1 
    % last file first
    thisout = fullfile(savedir, [srcf(end).name(1:end-4),fend]);
    thatout = fullfile(savedir, [srcf(1).name(1:end-4),fend]);
    plx_splitOut(zipname, splits(end), thatout, thisout);
    fprintf('\t\t%s\n',thisout)
    fprintf('\t\t%s\n',thatout)

else
    % last file first
    thisout = fullfile(savedir, [srcf(end).name(1:end-4),fend]);
    plx_splitOut(zipname, splits(end), ['part',num2str(length(splits)),'_',zipname], thisout);
    fprintf('\t\t%s\n',thisout)

    % all the ones inbetween
    for i = fliplr(1:length(splits)-1)
        if i > 1
            % next file to lop off the end
            thisout = fullfile(savedir, [srcf(i+1).name(1:end-4),fend]);
            plx_splitOut(['part',num2str(i+1),'_',zipname], splits(i), ['part',num2str(i),'_',zipname], thisout);
            fprintf('\t\t%s\n',thisout)
            delete(['part',num2str(i+1),'_',zipname])
        else
            % save the final pairing as originals
            thisout = fullfile(savedir, [srcf(i+1).name(1:end-4),fend]);
            lastout = fullfile(savedir, [srcf(i).name(1:end-4),fend]);
            plx_splitOut(['part',num2str(i+1),'_',zipname], splits(i), lastout, thisout);
            
            fprintf('\t\t%s\n',thisout)
            fprintf('\t\t%s\n',lastout)
            % clean up the partial junk
            delete(['part',num2str(i+1),'_',zipname])
            
            
        end
    end
end

fprintf('\t\tDone.\n')

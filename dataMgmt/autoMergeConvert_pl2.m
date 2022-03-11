function autoMergeConvert_pl2(baseDir, loadDir, useGUI, contFilterMode)
% function autoMergeConvert_pl2(baseDir, loadDir, useGUI, contFilterMode)
% 
% Recommend running without any inputs:
% Set current working directory to folder containing a dir for each recording
% day, formatted as  ./yyyymmdd, with original spike data files in ./yyyymmdd/spk
% Will present gui file selection windows for choosing which day, which files to sort, and corresponding chanMap to use
% 
% INPUTS - baseDir -- base directory, e.g. /home/huklab/
%          loadDir -- input directory where the pl data is
%          useGUI  -- flag to disable GUIs for scripting (def=1, show GUI file pickers where necessary)
%          contFilterMode -- 'CMR' default & recommended (...outlier versions 1-2x slower)
%                     - 'CAR', 'carOutliers', 'CMR', 'cmrOutliers', 'none'
%                     - "outliers" versions fill outliers across temporal dim (spikes) with nans
%                       before computing mean or median across channels.
%                     - only removes outliers from stat, not from actual matrix values
%                     - no outlier removal for temporal stat; at 20-40kHz, difference insignificant for any reasonable chunk size
% 
% NOTE ON CHANMAP FILES:   Needs extended version of chanMap files
%   -- All regular Kilosort vars (chanMap, chanMap0ind...xcoords, ycoords)
%   -- Additional 'czuba' fields of:
%       [name]  Name of chanMap file
%       [fs]    sampling frequency in Hz (e.g. fs = 40000;)
%       [configFxn]  function handle to Kilosort2 config file that should be used
%                    (e.g.  configFxn = str2func('nameOfMyConfig_KiloSort2'); % must be in path )
% 
% see also:  runKiloSort2, plx2raw
% 
% 2019-xx-xx  TBC  Cobbled together from prior Huk Lab Kilosort1 functions
% 2019-08-20  TBC  Cleaned & commented. Improved config & chanMap selection
% 2022-01-19  TBC  Created from old autoSort parts
% 2022-01-31  TBC  Converted old [timeRange] input to contFilterMode (def='CMR')
% 


%% Parse the inputs
if nargin<1 || isempty(baseDir)
    if exist('~/MLtoolbox', 'dir')
        baseDir = '~/MLtoolbox';
    end
end

if nargin<2 || isempty(loadDir)
    loadDir = pwd;
end

if nargin<3 || isempty(useGUI)
    useGUI = 1;
end

contFilters = {'CAR', 'carOutliers', 'CMR', 'cmrOutliers','none'};
if nargin<4 || isempty(contFilterMode)
    % continuous data filter mode (def='CMR' common median referencing across channels)
    % Filter Modes:  'CAR', 'carOutliers', 'CMR', 'cmrOutliers'
    contFilterMode = 'CMR';
elseif ~ischar(contFilterMode) || (ischar(contFilterMode) && ~any(strcmpi(contFilterMode, contFilters)))
    error(sprintf('Continuous data filter mode [contFilterMode] not recognized\nMust be one of the following: \t\t%s,\t\t%s,\t\t%s,\t\t%s\n',contFilters{:}));
end

% timeRange restriction not supported (messy & ineffective)
timeRange = [];

ok = true;

%% add necessary paths
% % addToPathWithoutGit( fullfile( baseDir, 'visbox','tbcToolbox') );
% % addToPathWithoutGit( fullfile( baseDir, 'Matlab Offline Files SDK') );  % path to plexon offline sdks
% % addToPathWithoutGit( fullfile( baseDir, 'npy-matlab') );                % path to npy-matlab scripts
% % addToPathWithoutGit( fullfile( baseDir, 'ks25') );                      % path to kilosort folder
% % addToPathWithoutGit( fullfile( baseDir, 'kilosort_utils') );            % path to kilosort utilities folder


if ~(exist(fullfile(loadDir,'spk'),'dir'))  % && exist(fullfile(loadDir,'pds'),'dir'))
    % Find all the dirs matching a recording session day (format:  'yyyymmdd')
    theseDirs=dir(loadDir);
    theseDirs=theseDirs([theseDirs.isdir]);
    theseDirs=theseDirs(~cellfun(@isempty,(regexp({theseDirs.name}, '^\d{8}$'))));
    
    if ~isempty(theseDirs) % && useGUI  % always start with dir gui
        [theseDirs, ok] = chooseFile(theseDirs, 'Select experiment day(s) to process: ', 'multiple');
    end
else
    % if [pwd] or [loadDir] is already a recording day directory
    % partition into outputs consistent with expectations (...krufty)
    [loadDir, theseDirs.name] = fileparts(loadDir);
    ok = true;
end

if isempty(theseDirs) || ~ok
    fprintf(2, 'No session dirs found in %s, or user canceled selection\n', loadDir);
    return
end


%% auto convert files
for iDir=1:length(theseDirs)
    thisdate=theseDirs(iDir).name;
    fprintf('~~~\n~~~\nStarting to convert session %s (%i of %i)\n',thisdate, iDir, length(theseDirs));
    
    % Get the spike files **Continuous data files should all be in ./spk/raw, not just ./spk**
    thispath = fullfile(loadDir, thisdate, 'spk', 'raw');
    b = dir( fullfile(thispath, '*plxC.pl2') );
    if isempty(b)   % no?? check direct/misnamed .pl2
        b=dir( fullfile(thispath, '*.pl2') );
        if isempty(b)   % check first level ./spk dir
            thispath = fileparts(thispath); %strip off 'raw' subdir
            b=dir( fullfile(thispath, '*plxC.pl2') );
        end
    end
    
    if useGUI==1
        if isempty(b)
            [b, ok] = uigetfile({'plx','pl2'}, 'Select spike data files for batch spike sorting: ');
        else
            [b, ok] = chooseFile(b, 'Select files for batch spike sorting: ');
        end
    else
        disp({b.name}')
    end

    if isempty(b) || ~ok
        warning('No spike files found in, or user canceled selection\n\t...skipping this dir:\t%s\n', thispath);
        continue
    end
    
    
    %% Merge & convert all selected files in session
    % compile each plxInfo(>>pI) & rawInfo(>>rI) struct, then generate merged at end of loop
    theseFiles = {b.name};

    % default raw file name for merged files
    jnk = strfind(theseFiles, thisdate);
    tt = jnk{1}+7;  %8==length(thisdate);
    rawFileOut = theseFiles{1}(1:tt);
    for i = 1:length(jnk)
        rawFileOut = [rawFileOut, theseFiles{i}(jnk{i}+8)]; %#ok<*AGROW> 
    end
    rawFileOut = [rawFileOut, theseFiles{1}(tt+2:end)];
    % (if present,) strip extension from raw filename string
    [~,rawFileOut] = fileparts(rawFileOut);
    %
    rawFileOut = inputdlg(sprintf('Create name for merged files:%s\n\nMerged file name (no extension):', sprintf('\n\t%s',b.name)), 'Merged File Name', 1, {rawFileOut});

    if isempty(rawFileOut)
        % user canceled, abort merge for this directory
        continue
    else
        % unpack cell output of inputdlg
        rawFileOut = rawFileOut{1}; 

        % For each spike file
        for iPL=1:length(b)
            % [thisFile] is the current spike file
            % [thispath] is the path to current spike file (...for [clumsy] reasons, could be in either ./spk or ./spk/raw)
            % [spkDir] is path to this session ./spk
            % [rawDir] is path to this session ./spk/raw  (converted data & "_rawInfo" file will always go here)

            [~,thisFile] = fileparts(b(iPL).name);
            fprintf('~~~\tProcessing file: \t%s\n\t...to: \t%s\n', thisFile, rawFileOut);

            if ~endsWith(thispath, 'raw')
                spkDir = fullfile(thispath);
                rawDir = fullfile(thispath,  'raw');
            else
                rawDir = thispath;
                spkDir = fileparts(rawDir);
            end
            if ~exist(rawDir, 'dir')
                mkdir(rawDir);
            end

            % ------------------------------------------------------------

            %% Compile plexon file info, and load strobe event timestamps & values

            % General PLX info struct
            % - 2nd input is flag to not plot events on this pass
            % - plx2raw will plot them when converting data
            plxInfo = getPlxInfo(fullfile(thispath, b(iPL).name), 0);

            % ------------------------------------------------------------
            % Find prospective channel map files ( <device>##ch<specs>##k.mat )
            cm = dir(fullfile(baseDir, 'kilosort_utils', sprintf('*%d*%dk*.mat', plxInfo.nChannels, round(plxInfo.freq/1000))));

            % copy channel map file
            if ~isempty(cm)
                if length(cm)>1
                    cm = chooseFile(cm, 'Select (1) channel map file: ', 'single');
                elseif ~isempty(cm)
                    % use newest match
                    [~,ii] = max([cm.datenum]);
                    cm = cm(ii);
                end
                chanMap = load(fullfile(cm.folder, cm.name));
            else
                warning('No chanMap match found.')
                chanMap = struct;
                keyboard
            end


            % ------------------------------------------------------------
            % Filter Modes:  'CAR', 'carOutliers', 'CMR', 'cmrOutliers'
            % - 'CMR' recommended
            % - outliers versions fill outliers across temporal dim (spikes) with nans
            %   before computing mean or median across channels.
            % - only removes outliers from stat, not from actual matrix values
            % - no outlier removal for temporal stat; at 20-40kHz, difference insignificant for any reasonable chunk size
            % contFilterMode = 'CMR';  %'carOutliers'; %'CMR'; % 'CMR'; %'none');   % 'CMR');

            % Continuous data source:   'SPKC'
            chanMap.contSource = 'SPKC';
            if isfield(chanMap,'contSource')
                contSource = chanMap.contSource;
            else
                contSource = [];
            end

            % [timeRange] restriction NOT RECOMMENDED/abandoned
            % - overcomplicates association with original data & event syncs
            %forceNew = 1;   timeRange = [1250, inf];

            if ~isempty(timeRange)
                % [timeRange] restriction NOT RECOMMENDED/abandoned
                % - overcomplicates association with original data & event syncs
                warning([mfilename,':timeRangeUnwise'], sprintf('[timeRange] restriction NOT RECOMMENDED/abandoned\n\t- overcomplicates association with original data & event syncs'));
                % append time range to output filename, if defined
                rawFileOut = sprintf('%s_t%04.0f-%04.0f', rawFileOut, timeRange);
            end


            %forceNew = 1;
            appendDat = iPL>1; % only append if more than one file, else plx2raw will default to create new & discard contents if exists

            % Convert raw data
            %if (forceNew || useGUI==2) || ~exist( fullfile( rawDir, [rawFileOut '.dat']), 'file')
            % create .dat raw file
            fprintf('Creating raw dat file from Plexon spikes\n');
            [plxInfo, rawInfo] = plx2raw(thispath, b(iPL).name, rawDir, contFilterMode, contSource, timeRange, rawFileOut, appendDat);   %   %'none');   % 'CMR');    %
            rI(iPL) = rawInfo;
            pI(iPL) = plxInfo;

        end % of loop through selected data files
        

        %% Create merged plxInfo struct
        % use first file as baseline info struct
        plxInfo = pI(1);
        for i = 2:length(pI)
            % concatenate strobes, shifting timestamps by previous file stop time
            plxInfo.strobes = [plxInfo.strobes; [pI(i).strobes(:,1) + plxInfo.rstop(end), pI(i).strobes(:,2)]];
            % rstart & rstop times
            plxInfo.rstart  = [plxInfo.rstart;  pI(i).rstart + plxInfo.rstop(end)];
            plxInfo.rstop   = [plxInfo.rstop;   pI(i).rstop + plxInfo.rstop(end)];
            % cumsum continuous sample counts (...shouldn't be used, but for posterity)
            % (all other columns MUST be identical)
            plxInfo.adNumCountFreqGains(:,2) = plxInfo.adNumCountFreqGains(:,2) + pI(i).adNumCountFreqGains(:,2);
            
        end
        plxInfo.nSamplesTot = max(plxInfo.adNumCountFreqGains(:,2));
        mergeInfo.plxInfo = pI;
        mergeInfo.rawInfo = rI;
        mergeInfo.source = theseFiles;  % cellstr(rI.paths.source); ...gah, krufty concatenation; just use cellstr from above
        % sync raw & plx info structs with paths to all relevant files
        rawInfo.paths.source = mergeInfo.source;
        plxInfo.paths = rawInfo.paths;
        % Update "_rawInfo.mat" file with mergedInfo struct
        save( rawInfo.paths.rawInfo, 'plxInfo', 'rawInfo', 'mergeInfo')

        
        %% auto copy channelmaps & rawInfo
        try
            % rawdatName = fullfile(rawDir, [rawFileOut '.dat']);
            outDir     = fullfile(spkDir, 'KiloSort', rawFileOut);

            if ~exist(outDir,'dir')
                mkdir(outDir);
            end

            % propagate chanmap & rawInfo files to output directory
            % - ks25 will also transfer these when saving, but redundancy doesn't hurt here
            copyfile(fullfile(cm.folder, cm.name), fullfile(outDir,'chanMap.mat'));
            copyfile(rawInfo.paths.rawInfo, outDir);

            fprintf('Kilo Sort output dir:\n\t%s\n', outDir);
        catch
            warning([mfilename,':ksOutputDirFail'], 'NOTE: unable to fully initialize default kilosort output dir:\t%s\n', outDir);
        end

        % Just convert raw data & save _rawInfo.mat, abort autosort so Kilosort GUI can be used
        fprintf('~!~\tRaw dat file conversion complete. Use Kilosort GUI to process raw file:\n\t\t%s.dat\n', rawFileOut);
        
    end

    fprintf('\nDone with directory:  %s\n\t%d of %d remaining.\n~~~~~~~~\n', thispath, length(theseDirs)-iDir, length(theseDirs))
    
end % n dirs

end % main function



% % % % % % % % %
%% Sub-Functions
% % % % % % % % %


%% chooseFile
function [fd, ok] = chooseFile(fd, titl, selMode)
    % fd == file struct output from:  fd = dir(____);
    % titl == file selection window title/prompt
    % selmode = 'single' or 'multiple'[default]
    %
    if nargin<3 || isempty(selMode) || strcmpi(selMode,'multiple')
        selMode = 'multiple';
    else
        selMode = 'single';
    end
    % select the last/most recent one by default
    defSel = length(fd);
    
    winsz = [400, 60+12*length(fd)];
    winsz(winsz>400) = 600;
    [fdsel, ok] = listdlg('liststring',{fd.name}, 'selectionmode',selMode, 'initialValue',defSel, 'listsize',winsz, 'Name', titl);
    if ok && ~isempty(fdsel)
        fd = fd(fdsel);
    else
        disp('User canceled selection')
        fd = [];
    end
end %end chooseFile

    function fcell = msg2fileCell(msg, starter)        % process return string to synced file/dirs
        msg = split(msg);
        fcell = msg(startsWith(msg, starter));
    end

% % % %% addToPathWithoutGit.m
% % % function addToPathWithoutGit(dir, excludes, withSubdirs)
% % %     if nargin<3 || withSubdirs
% % %         a = genpath(dir);
% % %         withSubdirs = ' and subdirectories';
% % %     else
% % %         a = dir;
% % %         withSubdirs = [];
% % %     end
% % %     
% % %     if isempty(a)
% % %         fprintf('%s not found...attempting to continue\n', dir);
% % %     else
% % %         b=textscan(a,'%s','delimiter',':');
% % %         b=b{1};
% % %         b(~cellfun(@isempty,strfind(b,'.git')))=[];
% % %         b(~cellfun(@isempty,strfind(b,'.svn')))=[];
% % %         if nargin>1
% % %             if ~iscell(excludes), excludes = {excludes}; end
% % %             for i = 1:numel(excludes)
% % %                 if ~isempty(excludes{i})
% % %                     b(~cellfun(@isempty,strfind(b, excludes{i})))=[];
% % %                 end
% % %             end
% % %         end
% % %         addpath(b{:})
% % %         disp([dir, withSubdirs, ' added to the path']);
% % %     end
% % % end %end addToPathWithoutGit

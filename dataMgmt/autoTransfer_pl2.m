function autoTransfer_pl2(baseDir, loadDir, useGUI, timeRange)
% INPUTS - baseDir -- base directory, e.g. /home/huklab/
%          loadDir -- input directory where the pl data is
% 
% ~~~ USE WITH CARE ~~~
% This function written to perform [somewhat esoteric] transfer of PLDAPS recording session data files & directories
% from external drive (fast ssd) to local hard drives
% - Creates destination file structure   <spkLocal>/<subj>/<yyyymmdd>
%   - <spkLocal> is your local fast storage drive (nvme, >=1GB/s)
%     - best use: create a default mountpoint at /spkLocal that points to this drive/directory
%     - do so across machines & analysis paths are consistent & transferrable
%   - <subj> is *hardcoded* below to 'kipp'
%   - <yyyymmdd> is the standard PLDAPS session directory name (e.g. 20220101)
% ~~~ HERE BE DRAGONS ~~~
% - /spkLocal/<subj>/<yyyymmdd>/spk/raw is created as a symlink tunnel to spkBulk
%   - <spkBulk> is your secondary larger storage volume
%     - follow same mountpoint recs described for spkLocal
%   - files in ./spk/raw will can be accessed as if they are present on spkLocal, but actually exist on spkBulk
%   - anything already in ./spk/raw on the source will be transferred to the symlink tunnel too (i.e. symlink tunnel is non-destructive)
%   - continuous data files from the ./spk direcotry will also be transferred into the ./spk/raw symlink tunnel
%     - detected based on all files matching:   "*plxC.pl2"   or   "*plxC.dat"
%
% - file transfers are performed with rsync
%   - /.phy directiories are ALWAYS EXCLUDED from sync
%     - just a bloated working directory of Phy; actual sorting/curated data & labels are stored as .npy files
% 
%
% Recommend running without any inputs. Set current working directory to folder containing a dir for each recording
% day, formatted as  ./yyyymmdd, with original spike data files in ./yyyymmdd/spk
% Will present gui file selection windows for choosing which day, which files to sort, and corresponding chanMap to use
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

if nargin<4 || isempty(timeRange)
    % use all available data by default
    timeRange = [];
end

ok = true;

%% Target directories
% subject dir
subj = 'kipp';
% Establish local copies of primary session files
spkLocal    = '/spkLocal';
% symlink continuous data files to bulk storage
% - separate drive, ideally large nvme, but could be HDD
spkBulk     = '/spkBulk';


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


%% Sync files from external drive to spkLocal & spkBulk
if ~startsWith(loadDir, spkLocal)

    for iDir=1:length(theseDirs)
        thisSession = theseDirs(iDir).name;
        fprintf('~~~\n~~~\tTransferring files from session %s (%i of %i)\n',thisSession, iDir, length(theseDirs));

        % create path component string for /<subj>/<yyyymmdd>/spk
        % - example:  fullfile( <spkLocal or spkBulk>, spkPath, <filename>);
        spkPath = fullfile(subj, thisSession, 'spk');

        %% Symlink ./spk/raw session directory to bulk storage drive
        % Initialize [symlinked] destination directories first
        if ~exist(fullfile(spkBulk, spkPath, 'raw'), 'dir')
            mkdir(fullfile(spkBulk, spkPath, 'raw'));
        end
        if ~exist(fullfile(spkLocal, spkPath), 'dir')
            mkdir(fullfile(spkLocal, spkPath));
        end
        rawPath = fullfile(spkLocal, spkPath, 'raw');
        if ~exist(rawPath, 'dir') || ~isSymlink(rawPath)
            % create symlink
            [err, msg] = system( sprintf('ln -sv %s %s', fullfile(spkBulk, spkPath, 'raw'), fullfile(spkLocal, spkPath )) );
            if ~err
                disp(msg)
            else
                fprintf(2, msg)
                keyboard
            end
        elseif isSymlink(rawPath)
            % do nothing, already setup to pass through
        else
            % unexpected scenario
            keyboard
            fprintf(2, '\n!!!\t ./spk/raw dir already exists in spkLocal:\n\t\t%s\n', rawPath)
            [~,statMsg] = system(sprintf('stat %s', rawPath));
            fprintf(2, '\n%s\n!!!\t PROCEED WITH CAUTION !!!\n', statMsg);
            keyboard

        end


        %% Copy NON-Continuous data (./eye, ./pds, ./spk >>exclude"*plxC.pl2" & "plxC.dat")
        % use system rsync to copy non-continuous files & directories to LOCAL storage drive
        [err, out] = system( sprintf(['rsync --exclude ".phy" --exclude "._*" '...          %standard rsync excludes
                                      '--exclude "*plxC.pl2" --exclude "*plxC.dat" '...     % exclude continuous data
                                      '--out-format %%n -Krluvhpt %s %s'], fullfile(loadDir, thisSession), fullfile(spkLocal, subj)) );   % send to session dir & subdirs to <spkLocal>/subj/  ("K" option follows raw symlink to spkBulk)

        % error check & process returned string to cellstr of synced file/dir names
        if ~err && ~isempty(out)
            disp(out)
            %             out = msg2fileCell(out, thisSession);
            %             fprintf('Synced the following files from %s to %s:\t(excl. continuous data)\n%s', loadDir, spkLocal, sprintf('\t%s\n',out{:}))
        else
            fprintf(2, out)
            keyboard
        end


        %% Copy Continuous data to spkLocal raw dir (which will pass through to spkBulk on its own)
        % FIRST, create ./spk/raw directory in [spkBulk]
        %   - create symlink w/in spkLocal pointing to spkBulk raw dir ("bulk data portal")
        % THEN, use system rsync to copy *continuous* data files & directories to BULK storage drive via the symlink in spkLocal
        contFilesSrc = fullfile(loadDir, thisSession, 'spk/');
        % Two-pass, non recursive
        % - ensure all continuous files from BOTH ./spk AND ./spk/raw subdirectory are compiled into ./spk/raw
        fprintf('\n~~~\tTransferring continuous data files from  %s  -->  %s \n', contFilesSrc, rawPath);
        [err, syncedC] = system( sprintf(['rsync --exclude ".phy" --exclude "._*" '... % --info=progress2 doesn't work with system() since msg is only returned on completion
                                          '--include "*plxC.pl2" --include "*.dat" --exclude "*" '...
                                          '--out-format %%n -rluvhpt %s %s'], contFilesSrc, rawPath) ); %*NOTE* trailing filesep tells rsync to only send files, not dir structure
        disp(syncedC);  if err, keyboard,   end
        if exist(fullfile(contFilesSrc,'raw'), 'dir')
            fprintf('\n~~~\tTransferring continuous data files from  %s  -->  %s \n', fullfile(contFilesSrc,'raw'), rawPath);
            [err, syncedC] = system( sprintf(['rsync --exclude ".phy" --exclude "._*" '... % --info=progress2 doesn't work with system() since msg is only returned on completion
                '--include "*plxC.pl2" --include "*.dat" --exclude "*" '...
                '--out-format %%n -rluvhpt %s %s'], fullfile(contFilesSrc,'raw/'), rawPath) ); %*NOTE* trailing filesep tells rsync to only send files, not dir structure
            disp(syncedC);  if err, keyboard,   end
        end

    end %for each session
end %transfer data



%% CONVERSION (broken out into separate file:  autoConvert_pl2)
% % % % update load dir to point to spkLocal for further data conversions
% % % origLoadDir = loadDir; %#ok<NASGU> % just in case..
% % % loadDir = fullfile(spkLocal, subj);
% % % 
% % % %% auto convert files
% % % for iDir=1:length(theseDirs)
% % %     thisdate=theseDirs(iDir).name;
% % %     fprintf('~~~\n~~~\nStarting to convert session %s (%i of %i)\n',thisdate, iDir, length(theseDirs));
% % %     
% % %     % Get the spike files
% % %     thispath = fullfile(loadDir, thisdate, 'spk', 'raw');
% % %     b = dir( fullfile(thispath, '*plxC.pl2') );
% % %     if isempty(b)   % no?? check direct/misnamed .pl2
% % %         b=dir( fullfile(thispath, '*.pl2') );
% % %         if isempty(b)   % check first level ./spk dir
% % %             thispath = fileparts(thispath); %strip off 'raw' subdir
% % %             b=dir( fullfile(thispath, '*plxC.pl2') );
% % %         end
% % %     end
% % %     
% % %     if useGUI==1
% % %         if isempty(b)
% % %             [b, ok] = uigetfile({'plx','pl2'}, 'Select spike data files for batch spike sorting: ');
% % %         else
% % %             [b, ok] = chooseFile(b, 'Select files for batch spike sorting: ');
% % %         end
% % %     else
% % %         disp({b.name}')
% % %     end
% % % 
% % %     if isempty(b) || ~ok
% % %         warning('No spike files found in %s, or user canceled selection\n\tWill attempt to continue anyway...\n', thispath);
% % %     end
% % %     
% % %     
% % %     
% % %     % For each spike file
% % %     for iPL=1:length(b)
% % %         % [thisFile] is the current spike file
% % %         % [thispath] is the path to current spike file (...for [clumsy] reasons, could be in either ./spk or ./spk/raw)
% % %         % [spkDir] is path to this session ./spk
% % %         % [rawDir] is path to this session ./spk/raw  (converted data & "_rawInfo" file will always go here)
% % % 
% % %         [~,thisFile] = fileparts(b(iPL).name);
% % %         fprintf('~~~\tProcessing file:  %s\n', thisFile);
% % %         
% % %         if ~endsWith(thispath, 'raw')
% % %             spkDir = fullfile(thispath);
% % %             rawDir = fullfile(thispath,  'raw');
% % %         else
% % %             rawDir = thispath;
% % %             spkDir = fileparts(rawDir);
% % %         end
% % %         if ~exist(rawDir, 'dir')
% % %             mkdir(rawDir);
% % %         end
% % %         
% % %         % (if present,) strip extension from raw filename string
% % %         [~,rawFileOut] = fileparts(thisFile);
% % %         % ------------------------------------------------------------
% % %         
% % %         %% Compile plexon file info, and load strobe event timestamps & values
% % %         
% % %         % General PLX info struct
% % %         % - 2nd input is flag to not plot events on this pass
% % %         % - plx2raw will plot them when converting data
% % %         plxInfo = getPlxInfo(fullfile(thispath, b(iPL).name), 0);
% % % 
% % %         % ------------------------------------------------------------
% % %         % Find prospective channel map files ( <device>##ch<specs>##k.mat )
% % %         cm = dir(fullfile(baseDir, 'kilosort_utils', sprintf('*%d*%dk*.mat', plxInfo.nChannels, round(plxInfo.freq/1000))));
% % % 
% % %         % copy channel map file
% % %         if ~isempty(cm)
% % %             if length(cm)>1
% % %                 cm = chooseFile(cm, 'Select (1) channel map file: ', 'single');
% % %             elseif ~isempty(cm)
% % %                 % use newest match
% % %                 [~,ii] = max([cm.datenum]);
% % %                 cm = cm(ii);
% % %             end
% % %             chanMap = load(fullfile(cm.folder, cm.name));
% % %         else
% % %             warning('No chanMap match found.')
% % %             chanMap = struct;
% % %             keyboard
% % %         end
% % %         
% % % 
% % %         % ------------------------------------------------------------
% % %         % Filter Modes:  'CAR', 'carOutliers', 'CMR', 'cmrOutliers'
% % %         % - 'CMR' recommended
% % %         % - outliers versions fill outliers across temporal dim (spikes) with nans
% % %         %   before computing mean or median across channels.
% % %         % - only removes outliers from stat, not from actual matrix values
% % %         % - no outlier removal for temporal stat; at 20-40kHz, difference insignificant for any reasonable chunk size
% % %         % contFilterMode = 'carOutliers'; %'CMR'; % 'CMR'; %'none');   % 'CMR');
% % %         contFilterMode = 'CMR';
% % %         
% % %         % Continuous data source:   'SPKC'  
% % %         chanMap.contSource = 'SPKC';
% % %         if isfield(chanMap,'contSource')
% % %             contSource = chanMap.contSource;
% % %         else
% % %             contSource = [];
% % %         end
% % %         
% % %         % [timeRange] restriction NOT RECOMMENDED/abandoned
% % %         % - overcomplicates association with original data & event syncs
% % %         %forceNew = 1;   timeRange = [1250, inf];
% % %         
% % %         if ~isempty(timeRange)
% % %             % append time range to output filename, if defined
% % %             rawFileOut = sprintf('%s_t%04.0f-%04.0f', rawFileOut, timeRange);
% % %         end
% % %         
% % %         
% % %         forceNew = 1;
% % %         
% % %         
% % %         % Convert raw data
% % %         if (forceNew || useGUI==2) || ~exist( fullfile( rawDir, [rawFileOut '.dat']), 'file')
% % %             % create .dat raw file
% % %             fprintf('Creating raw dat file from Plexon spikes\n');
% % %             rawInfo = plx2raw(thispath, b(iPL).name, rawDir, contFilterMode, contSource, timeRange, rawFileOut);    %#ok<NASGU>    %'none');   % 'CMR');    %
% % %             
% % %         else
% % %             rawInfo = fullfile( rawDir, [rawFileOut,'_rawInfo.mat']);
% % %             if exist( rawInfo, 'file')
% % %                 fprintf('Found existing raw dat file; loading rawInfo struct...\n');
% % %                 rawInfo = load(rawInfo); %#ok<NASGU> 
% % %                 
% % %             else
% % %                 warning('Raw file for %s already exists, but no info struct was saved for it, so I don''t know whats in there.\n Raw file should be recreated...', thisFile)
% % %                 keyboard
% % %                 
% % %             end
% % %         end
% % %         
% % %         %% auto copy channelmaps
% % %         % rawdatName = fullfile(rawDir, [rawFileOut '.dat']);
% % %         outDir     = fullfile(spkDir, 'KiloSort', thisFile);
% % %         
% % %         if ~exist(outDir,'dir')
% % %             mkdir(outDir);
% % %         end
% % %         
% % %         copyfile(fullfile(cm.folder, cm.name), fullfile(outDir,'chanMap.mat'));
% % %         
% % %         fprintf('Kilo Sort output dir:\n\t%s\n', outDir);
% % % 
% % %         
% % %         %fprintf('~!~\tRaw dat file conversion complete. Use Kilosort GUI to process raw file:\n\t%s.dat\n', rawFileOut);
% % %         %return
% % %         
% % %         
% % %         
% % %         if 1 %useGUI==2
% % %             % Just convert raw data & save _rawInfo.mat, abort autosort so Kilosort GUI can be used
% % %             fprintf('~!~\tRaw dat file conversion complete. Use Kilosort GUI to process raw file:\n\t%s.dat\n', rawFileOut);
% % %             
% % %             % % %         else
% % %             % % %             %% Actual call to KiloSort!
% % %             % % %             % ----------------------------
% % %             % % %             % This code is long out of date...needs update for "ks25" version of Kilosort
% % %             % % %             % ----------------------------
% % %             % % %             if ~exist(fullfile(outDir, 'spike_times.npy'), 'file')
% % %             % % %                 rez = runKiloSort2(outDir, thisFile, 'pl2');
% % %             % % %                 fprintf('\n\tProcessing time: %3.1f min\n\t%s done.\n', toc/60, thisFile);
% % %             % % %             else
% % %             % % %                 fprintf(2, '\tKilosort has already been run for this file/path. Moving on.\n');
% % %             % % %             end
% % %             % % %
% % %             % % %             % Make a guess at labeling output figures....
% % %             % % %             H = gcf;
% % %             % % %             try                             %#ok<*TRYNC>
% % %             % % %                 set(gcf, 'name',[thisFile,'_templates'], 'tag',fullfile(rez.ops.saveDir,'figs'));
% % %             % % %                 figureFS(gcf, 'landscape',[8,5]);
% % %             % % %             end
% % %             % % %             try
% % %             % % %                 figure(H.Number-1);
% % %             % % %                 set(gcf, 'name',[thisFile,'_batch2batch'], 'tag',fullfile(rez.ops.saveDir,'figs'));
% % %             % % %                 figureFS(gcf, 'landscape',[8,5]);
% % %             % % %             end
% % %             % % %
% % %         end % ksGUI breakout        
% % % 
% % %     end % n pl files within dir
% % %     fprintf('\nDone with directory:  %s\n\t%d of %d remaining.\n~~~~~~~~\n', thispath, length(theseDirs)-iDir, length(theseDirs))
% % %     
% % % end % n dirs

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

function dv = initDaily(baseName, basePath, spkSrc, unitArgs)
% Initialize structure & file paths for daily recording files with--unit ordered--cell array style spike structures
% and parallel toolbox for loading of spike data.
%   (based on code infrastructure written by T.Czuba at Einstein, 2012-2016)
% 
%       Array data and thomas tetrode data now referred to as "arr" and "tom"
%           (instead of "nev" and "plx" since file type is no longer exclusive)
%       ? Should a spkSrc distinction be made for uprobe data recorded on MAP system vs Omniplex?
% 
%       Data structure formatted as:
%           id:     unit id (chan*10 + unit#)
%           ts:     timestamp cell (nUnits x 1; in msec)
%           sync:   sync times  (in msec)
%           snr:    waveform snr (nUnits x 1)
%           wf:     struct with mean and quantiles of waveforms (fields mu & ci)
%                   * can also request to have .all field with all waveforms...will make file huge & not recommended
%           n:      nSpikes per unit (nUnits x 1)
% 
% ....updating in progress Sept. 2014. Working with just array data [arr], so code to get data from thomas system [tom]
%     (e.g. tetrodes) is not yet updated.

%% Parse inputs

% Which spikes to load/use
if nargin<4 || isempty(unitArgs)
    unitArgs = [666, 0]; % Default:  666=load all spikes from PLX/PL2 file, collapse all sorted & unsorted
end

if isnumeric(unitArgs) && length(unitArgs)<2
    unitArgs = [unitArgs, 0];   % don't get wfs by default

elseif isstring(unitArgs) || ischar(unitArgs)
    unitArgs = string(unitArgs); % ensure string class
    % For kilosort outputs, kilosort dir should be in /<sessionDate>/spk/KiloSort/<kilosortOutputDir>
    %   unitArgs = "kilo";  will find kilosort output dir that matches plexon file name
    %   unitArgs = "<name of kilosort output dir>"
end


%% Identify Spike Source
if nargin<3 || isempty(spkSrc)
    spkSrc = {'uprb'};
elseif ischar(spkSrc) || iscell(spkSrc)
    spkSrc = cellstr(spkSrc);
else
    error(fprintf('ERROR:\tWrong input type for spike source(s);  ''spkSrc''.\n\tShould be string (or cellstr) of spike source abreviations to be initialized.\n\t(abrev. should match directory names & follow convention; e.g. "arr", "tom", "uprb"...etc)\n'))
end


% Base Path (usually the recording day directory)
if nargin<2 || isempty(basePath)
    basePath = pwd;
    basePath = homeTilda(basePath);
end
dv.paths.root = basePath;

% PDS file selection modal if [baseName] not specified
if nargin<1 || isempty(baseName)
    fd = dir(fullfile(dv.paths.root, 'pds', '*.PDS'));
    fd = chooseFile(fd, 'Select PLDAPS data file:');
    % split name into parts
    fdParts = strsplit(fd.name, '_');
    % make baseName
    baseName = sprintf('%s_%s', fdParts{1:2});
end

%% Get PLDAPS data structure
pdsPath = fullfile(basePath, 'pds');
% check for pre-compiled PDS struct first
fd = dir( fullfile(pdsPath, ['*',baseName '*.mat']) );
if isempty(fd)
    fd = dir( fullfile(pdsPath, ['*',baseName '*.PDS']) );
end

% Load PDS file
if ~isempty(fd)
    if length(fd)>1
        fd = chooseFile(fd, 'Select PDS source file:');
    end
    
    switch lower(fd.name(end-2:end))
        case 'mat'
            % use pre-compiled data struct
            thisfile = fullfile(pdsPath, [fd.name(1:end-3),'mat']);
            fprintf('Loading pre-compiled analysis struct:\n\t%s\n', thisfile);
            dv = load( thisfile );
            
        case 'pds'
            % load PDS struct into dv output struct
            dv.paths.pds = fullfile(pdsPath, fd.name);
            dv.pds = pdsImport(dv.paths.pds);
    end
    
    % PDS stimulus file info
    [~, info.pdsName] = fileparts(dv.pds.baseParams.session.file);
    info.stimType = {dv.pds.condMatrix.modNames.currentStim{1}, dv.pds.baseParams.session.caller.name};
    info.viewDist = dv.pds.baseParams.display.viewdist;
    dv.info = info;
end


%% Find spike file(s)
% Spikes are ALWAYS loaded in different function. This just searches for proper files
% and establishes paths

% Start with expected naming convention
seshTime = datestr(dv.pds.baseParams.session.initTime, 'yyyymmdd');
if isstring(dv.pds.baseParams.session.subject)
    subjIni = dv.pds.baseParams.session.subject(1);
    subjDescriptor = char(dv.pds.baseParams.session.subject(2));
    spkBase = sprintf('%s_%s%s', subjIni, seshTime, subjDescriptor(1));
else
    % subject input was simple '' char, not "" string class
    subjIni = dv.pds.baseParams.session.subject;
    subjDescriptor = [];
    spkBase = sprintf('%s_%s%s', subjIni, seshTime, subjDescriptor);
end

% Warn if likely to fail finding a matching spike file
if ~contains(dv.paths.pds, spkBase)
    fprintf(2, '\nPossible mismatch between the spike file base:\t\t%s\nand PDS file name:\t\t%s\n', spkBase, dv.paths.pds)
    fprintf('File-letter may not have been properly incremented when PLDAPS file was run.\nCheck first character of <subjDescriptor> variable in %s.m...try:\n', mfilename)
    fprintf('\n\tsubjDescriptor = [char(subjDescriptor(1)*1+1), subjDescriptor(2:end)];\n\tspkBase = sprintf(''%%s_%%s%%s'', subjIni, seshTime, subjDescriptor(1));\n\n')
    % There's a helper function (fixPdsFilename.m) to update PDS file name & relevant internal struct components,
    % but post-hoc changes to source data files is dangerous and not recommended. ...do it right the first time.
    keyboard
end

%% For each spike source...
for ss = 1:length(spkSrc)
    
    % ...currently assuming only a single "spk" dir exists with spike data for this file
    % thispath = fullfile(basePath, spkSrc{ss}, 'plx');
    % thispath = fullfile(basePath, spkSrc{ss}, 'sort');
    thispath = fullfile(basePath, 'spk');

    if isnumeric(unitArgs) && unitArgs(1)<0
        % force ignore sorted results
        fd = dir( fullfile(thispath,['*',spkBase,'*plx.pl2']) );
        unitArgs(1) = abs(unitArgs(1));
    else
        % first try to get sorted file (Offline Sorter)
        fd = dir( fullfile(thispath, 'sort', ['*',spkBase,'*sort*_plxC.pl2']));
        if ~isempty(fd)
            % use most recent sort file
            fd = fd(end);
            thispath = fd.folder;
            fprintf('~~~\tFound sorted %s file(s):\n\t\t%s\n', spkSrc{ss}, fullfile(thispath, fd.name));
        end

    end
    
    % if that found nothing, check for raw data
    if isempty(fd)        
        % Prioritize .various plx file types (...a moving target)
        fd = dir( fullfile(thispath,['*',spkBase,'*plx.pl2']) ); % initial plx file converted to pl2 (MUCH faster to load & handle)
        if isempty(fd)
            fd = dir( fullfile(thispath,['*',spkBase,'*.pl2']) ); % any other matching .pl2 files?
            if isempty(fd)
                fd = dir( fullfile(thispath,['*',spkBase,'*.plx']) ); % regular old plx file
                if isempty(fd)
                    fd = dir( fullfile(thispath,['*',baseName,'*.pl2']) ); % try the input baseName                
                    if isempty(fd)
                        fd = dir( fullfile(thispath,['*',baseName,'*.plx']) ); % regular old plx file matching baseName??....last chance!!
                    end
                end
            end
        end
        if isempty(fd)
            fprintf(2, '***\tCouldn''t find likely matches for %s file:\n\t\t%s\n***\t...attempting to continue without it.\n\n', spkSrc{ss}, fullfile(thispath, ['*',baseName,'*'])); %,'/raw/'
        else
            fprintf('Found %s spike data file(s):\n\t%s\n', spkSrc{ss}, fullfile(thispath, fd.name));
        end
    end
    
    % Assign file
    if ~isempty(fd)
        if length(fd)>1
            fd = chooseFile(fd, sprintf('Select %s source file:',spkSrc{ss}));
        end
        dv.paths.(spkSrc{ss}) = fullfile(thispath, fd.name);
        dv.(spkSrc{ss}) = struct();
    end

    
    %% Create PreProcessed name & load if exists
    % file name for saving synced spike output (makes preprocessed subset of jumbo data files for easier loading)
    [pdsPath, pdsFile, pdsExt] = fileparts(dv.paths.pds);
    [spkPath, spkFile, spkExt] = fileparts(dv.paths.(spkSrc{ss}));
    % Subdir for preprocessed data
    if ~exist(fullfile(spkPath, 'sync'), 'dir')
        mkdir(fullfile(spkPath,'sync'));
    end
    
    if isnumeric(unitArgs)
        modeStr = sprintf('mode%d', unitArgs(1));
    else
        modeStr = sprintf('mode-%s', unitArgs(1));
    end
    
    syncedOutput = fullfile(spkPath, 'sync', sprintf('%s_%s_%s.mat', spkFile, pdsFile, modeStr));
    
    if exist(syncedOutput, 'file')
        dv.(spkSrc{ss}) = load(syncedOutput);
        fprintf(2, '~~~\tPre-processed & synced <%s> output struct loaded from:\n\t\t%s\n', spkSrc{ss}, syncedOutput);
        
    else
        %% Sync spike & PDS clocks (only loads event timestamps, so should always be quick)
        % Synchronize PDS & Plexon file clocks and trim to valid trials
        switch lower(dv.paths.(spkSrc{ss})(end-2:end))
            case {'plx','pl2'}
                % This sync applies to all PLDAPS-Plexon based experiment sessions
                % -- It only operates on plexon event channels & PTB trial times, so outputs apply
                %    regardless of sort method (unsorted, Plexon Offline Sorter, KiloSort...etc.)
                sync = syncPlexon2PDS(dv.pds, dv.paths.(spkSrc{ss}) );
                
                
                if isnumeric(unitArgs)
                    % load spikes from plx/pl2 file (unsorted waveforms or Offline Sort outputs)
                    dv.(spkSrc{ss}) = plx_readerPar_Pldaps(dv.paths.(spkSrc{ss}), unitArgs, sync);
                    
                    if isempty(dv.(spkSrc{ss}).info.comment)
                        dv.(spkSrc{ss}).info.comment = input(sprintf('Info string for this file (e.g.  "[channels]:  [gridLoc]-[depth in microns]" )\n\t:: '), 's');
                    end
                    
                else
                    % Determine path to kilosort output directory
                    switch unitArgs(1)
                        case 'kilo'
                            % use kilosort output directory matching this plx file
                            kilopath = getKiloPath(dv.paths.(spkSrc{ss}));
                        otherwise
                            % use different kilosort output directory (e.g. merged outputs of multiple plx/pl2 files)
                            kilopath = getKiloPath(dv.paths.(spkSrc{ss}), unitArgs(1));
                    end
                    
                    % Generate the spk output struct (finally!)
                    dv.(spkSrc{ss}) = getSpikes_kilo(kilopath, sync);
                    
                    % TODO: Ensure any time conversions applied to syncs are also applied to spikes
                    %       (see plx_readerPar_Pldaps.m)

                    
                end
                               
            otherwise
                dv.(spkSrc{ss}).sync = [];
        end
        
        
        dv.(spkSrc{ss}).info.preprocessed = syncedOutput;
        dv.(spkSrc{ss}).info.calcdate = datestr(now, 31);
        % save a pre-processed & synced mat file for faster loading in the future
        theseSpikes = dv.(spkSrc{ss});
        eval(sprintf('save(''%s'', ''-struct'',''theseSpikes'');', dv.(spkSrc{ss}).info.preprocessed));
        fprintf('~~~\tPre-processed & synced [.%s] spike output struct saved as:\n\t\t%s\n', spkSrc{ss}, theseSpikes.info.preprocessed );

    end
   
end


%% Eye tracking file path
try
    eyeFilename = [dv.pds.baseParams.session.file(1:end-3), 'edf'];
    if exist(fullfile(dv.paths.root, 'eye', eyeFilename), 'file')
        % relative path
        dv.paths.eye = fullfile('.', 'eye', eyeFilename);
    end
end

if ~evalin('caller', sprintf('exist(''figDir'',''var'') && ~isempty(figDir)'))
    figDir = fullfile(dv.paths.root, 'figs');   %'~/Dropbox/Science/projects/3dViewDist/kipp/figs';
    evalin('caller', sprintf('figDir = ''%s''', figDir))
end



%% DONE.
end

% % % % % % % % %
% Sub-functions
% % % % % % % % %

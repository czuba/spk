function filename = fixPdsFilename(filebase)

% file name
if nargin<1 || isempty(filebase)
    filebase = []%
%     defarg('filebase', sprintf('%s_%s', 'kipp', datestr(now,'yyyymmdd')))
%     filebase = cell2mat(inputdlg('rfPos File Name String     (e.g. ''<subj>_YYYYMMDDa_<gridLoc>'')', 'FileName', 1, {filebase}));
end


fd = dir(fullfile(pwd, ['*',filebase,'*.PDS']));
if isempty(fd)
    error('No PDS files detected in this directory.\n%s must be run from within the dir containing PDS files that need to be updated.',mfilename)

elseif numel(fd)>1
    fd = chooseFile(fd, 'Select PLDAPS data file:', 'single');
    %     % split name into parts
    %     fdParts = strsplit(fd.name, '_');
    %     % make baseName
    %     baseName = sprintf('%s_%s', fdParts{1:2});
end

origName = fd.name;
bkupName = sprintf('%s.bkup%s', fd.name, datestr(now,'YYYYMMDD'));
% create backup (this carries original 'date modified' over to the backup)
movefile(origName, bkupName);
copyfile(bkupName, origName);

% only load fields that determine file name
pds = load(origName, '-mat', 'baseParams');

fprintLineBreak;
fprintf('Backup PDS saved at:\n\t%s\n', bkupName);


%% Update filename params
curFileName = pds.baseParams.session.file;
sessionDateStr = datestr(pds.baseParams.session.initTime, 'yyyymmdd');
sessionTimeStr = datestr(pds.baseParams.session.initTime, 'HHMM');
subjStr = pds.baseParams.session.subject; % subject name
labelStr = pds.baseParams.session.experimentSetupFile; % session label

if isstring(subjStr)
    s0 = cellstr(subjStr);
    % poll user for corrections
    s = inputdlg(s0 , sprintf('Current file name:\t%s',curFileName), 1, s0);
    % Apply to original, matching input dimensions
    subjStr = reshape(string(s), size(pds.baseParams.session.subject));
    pds.baseParams.session.subject = subjStr;
else
    error('Handling for non string input not coded yet.')
    
end


%% Make file string
% if subject is string array (Matlab >2017), use second element as
if isa(subjStr, 'string')
    labelStr = [labelStr, char(subjStr(1,2))];
    subjStr = subjStr(1,1);
end
% ...clean up formatting for file name
if ~isempty(subjStr)
    subjStr = sprintf('%s_', subjStr);
end
if ~isempty(labelStr)
    labelStr = sprintf('%s', labelStr);
end

% Session filename
filename = sprintf('%s%s%s_%s.PDS',...
    subjStr,...
    sessionDateStr,...
    labelStr, ...
    sessionTimeStr);

%% apply filename string to original
pds.baseParams.session.file = filename;

% ~~~ NO!  Don't do this! ~~~
% pdsCore should remain an untouched record of running state
% % % % Apply to pds internals (PDS.pdsCore)
% % % pds.pdsCore.initialParameters{find(strcmpi(string(pds.pdsCore.initialParameterNames), 'SessionParameters'))}.session.file % = filename
% ~~~

%% Report changes & return
fprintf('Updates saved to:\n\t%s\n', origName)
fprintf(2, 'Consider manual updating to:\n\t%s\n', filename);
fprintLineBreak;

save(origName, '-append', '-struct','pds')
end


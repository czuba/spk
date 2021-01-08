function [evStruct] = getPlxEvents(plxPath)
% function [strobes] = getPlxEvents(plxFilename)
% 
% Loads all strobed event syncs from [plxPath], results are returned
% in fields of evStruct.
% 
% If [plxPath] input is a plxInfo struct, event output fields will be
% appended & all will be returned collectively as [evStruct]
% 
%
% INPUTS:
%   plxPath     = Full path to Plexon data file (.plx or .pl2) 
% 
% OUTPUTS:
%   evStruct fields:
%   .path           = full path to source file
%   .strobes        = n-by-2 vector of [timestamp, value]
%   .evNames        = event channel names (as cellstr; ommit any empty evChans?)
%   .rstart         = timestamps of rstart events
%   .rstop          = timestamps of rstop events
% 
%   timestamp conversion functions:
%   - relative to data file creation time
%   - For qualitative use only; plx file start time only precise to seconds (?!?)
%   .plx2clock      = convert to clock time; * outputs "duration" class for readability
%                   e.g. evClock = evStruct.plx2clock(evStruct.strobes(:,1));
% 
% Note:
%   Timestamp values are returned just as they are in the plx file, thus may not immediately
%   coincide with values returned from syncPlexon2PDS.m (or other PLDAPS analysis functions)
%   that are typically synced to stimulus (psychtoolbox) computer timestamps.
%   
% 
% see also:  getPlxInfo, plx_readerPar_Pldaps, syncPlexon2PDS
% ---
% 2020-01-06  TBC  Wrote it.
% 

%% Parse inputs
evStruct = struct;  %('path',[],'strobes',[], 'evNames',[],'rstart',[],'rstop',[],'plx2clock',[]); % init output struct

if nargin <1 || isempty(plxPath)
    [fileName, filePath] = uigetfile('*.*', 'Select spike source file');
    if ~fileName
        % User pressed cancel
        evStruct = [];
        return
    end
    plxPath = fullfile(filePath, fileName);
    
elseif isstruct(plxPath)
    % input was plxInfo struct
    evStruct = plxPath;
    plxPath = plxPath.path;
end

% expand home directory path (else will crash in mexplex [Plexon SDK])
if strcmp(plxPath(1),'~')
    [~,jnk] = system('echo $HOME');  jnk = jnk(1:end-1);
    plxPath = fullfile(jnk,plxPath(2:end));
end

if ~isfield(evStruct,'path')
    evStruct.path = plxPath;
end

%% Load events
% Get strobed plexon events (Omniplex uses 'channel 257' for the strobbed events)
[~, plxStrobeTs, plxStrobeVals] = plx_event_ts(plxPath, plx_event_resolve_channel(plxPath, 'Strobed'));

if plxStrobeTs == -1
    fprintf(['\tPlexon MAP system detected. Will use hardcoded channel #257 for strobed sync signals.\n'])
    [~, plxStrobeTs, plxStrobeVals] = plx_event_ts(plxPath, 257);
end
fprintf('\t%2.2f sec to load events from plx file:\t%s\n\n', toc, plxPath);

evStruct.strobes = [plxStrobeTs, plxStrobeVals];

%% Post-process & info
% get names of each event channel
[~, evNames] = plx_event_names(plxPath);
evStruct.evNames = cellstr(evNames);

% RSTART & RSTOP
try
    % This will only work with PL2 files
    % ...typical PlexonSDK functions are borked for RSTART & RSTOP events (circa 2020)
    ev = PL2EventTs(plxPath, 'RSTART');
    evStruct.rstart = ev.Ts;
    ev = PL2EventTs(plxPath, 'RSTOP');
    evStruct.rstop = ev.Ts;
catch
    fprintf(2,'!~!\tWarning: unable to read RSTART & RSTOP events from plx file\n')
    evStruct.rstart = [];
    evStruct.rstop = [];
end

% Timestamp to clock conversion
if isfield(evStruct, 'dateStr')
    % use value(s) from input plxInfo struct
    t0 = evStruct.dateStr;
else
    % get dateStr directly
      [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, Duration, t0] = plx_information(plxPath);
end

% convert to seconds of day
t0 = datetime(t0, 'InputFormat','MM/dd/uuuu HH:mm:ss');
t0time = timeofday(t0);
evStruct.plx2clock = @(x) duration(t0time + seconds(x), 'format','hh:mm:ss.SS');

% Special case: PLDAPS syncs detected
if all(evStruct.strobes(1:6,2)>2^14)
    fprintf('\tApplying [qualitative] sync of event times to PLDAPS clock\n')
    % - First set of 6 sync values sent by PLDAPS are "unique trial identifiers":
    %   [trial#, month, day, hour, minute, second]
    % - values are shifted into the upper part of strobed word register by adding 2^14
    % - use this time stamp as [coarse] sync to stimulus computer
    % - - !! Note: this is NOT a substitute for proper data syncing (via syncPlexon2PDS.m) !!
    t0prime = datetime([year(t0), evStruct.strobes(2:6,2)'-2^14]); % strobed datetime of first sync(s)
    s1 = evStruct.plx2clock(evStruct.strobes(1,1)); % plx datetime of first sync
    % update t0time & conversion function handle
    t0time = t0time + (timeofday(t0prime)-s1);
    evStruct.plx2clock = @(x) duration(t0time + seconds(x), 'format','hh:mm:ss.SS');%x + t0sec;
end

% clock time of data acquisition start [approx.]
evStruct.t0time = t0time;

end %main function


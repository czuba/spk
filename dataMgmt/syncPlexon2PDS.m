function [nsync, PDS] = syncPlexon2PDS(PDS, plxFilename, verbose)
% function nsync = syncPlexon2PDS(PDS, plxFilename, verbose)
% Sync Psychtoolbox, Datapixx, & Plexon clocks
% INPUTS:
%   [PDS]   could be:
%           PDS filename in standard location: <pwd>/pds/<filename>
%           Fullpath to PDS file
%           or preloaded PDS struct
%   [plxFilename] if not input, will attempt to find based on PDS filename
%   [verbose]   Flag to create clock syncing plots (optional, default=1)
% 
% OUTPUT: [nsync] struct (see end of fxn for complete output struct contents)
%                 Primary fields: .strobe, .plxTrialTs, .info
%                 Conversion fxns: .plx2ptb, .ptb2plx, .plx2dpx, .dpx2plx
% 
%         [PDS] data struct loaded with pdsImport.m
%         -- adds basic info fields to pds data struct (currently Uprobe centric)
% 
% Use conversion fxn handles to convert timestamps
% e.g.  spikeTs = nsync.plx2ptb(spk.ts{thisUnit})
% 
% ***Must have Plexon SDK Matlab files are in your path***
% 
% 
% Originally: jk wrote the underpinnings, eh started commenting,
%             tc further/partially decyphered....
% 
% 2018-05-29  TBC   Somewhat scrubbed/commented version of sync code used with old 8-bit-centric
%                   version of PLDAPS event syncing.
% 2019-08-16  TBC   This got a major overhaul during spring VSS prep.
% 2019-08-16  TBC   Cleaned up a bunch of detritus...see one git version back to see commented out code
% 

%% Parse inputs
if nargin < 3
    verbose = 1;
end

% if PDS input is filename, load PDS struct
if ~isstruct(PDS) 
    if exist(PDS, 'file')
        pdsFilename = PDS;
    elseif exist(fullfile('.','pds',PDS), 'file')
        pdsFilename = fullfile('.','pds',PDS);
    end
    tic
    PDS = pdsImport(pdsFilename);
    fprintf('\t%2.2f sec to load PDS file:\t%s\n', toc, pdsFilename);

else
    try
        pdsFilename = PDS.baseParams.session.file;
    catch
        % old wordy style PDS file
        pdsFilename = PDS.initialParametersMerged.session.file;
    end
end

[basePath, jnk] = fileparts(fileparts(pdsFilename));
if ~strcmpi(jnk,'pds')
    % detect non-standard session structure and try to cope
    basePath = fullfile(basePath, jnk);
end

if nargin < 2
    % assume PDS & Plexon file should have corresponding filenames
    try
        sessionDateStr = datestr(PDS.baseParams.session.initTime, 'yyyymmdd');
        if isstring(PDS.baseParams.session.subject)
            subjStr = PDS.baseParams.session.subject(1);
            if length(PDS.baseParams.session.subject)>1
                labelStr = PDS.baseParams.session.subject(2);
            end
        else
            subjStr = PDS.baseParams.session.subject;
            labelStr = p.defaultParameters.session.experimentSetupFile;
        end
        plxFilename = fullfile(basePath,'spk', sprintf('%s_%s%s.plx', subjStr,sessionDateStr,labelStr));
        
    catch
        % old wordy style PDS file...no guarantees
        subjStr = PDS.initialParametersMerged.session.subject;
        sessionDateStr = datestr(PDS.initialParametersMerged.session.initTime, 'yyyymmdd');
        labelStr = PDS.initialParametersMerged.session.experimentSetupFile;
        plxFilename = fullfile('.',sprintf('%s_%s%s.plx', subjStr,sessionDateStr, '*')); % trim final '_HHMM.PDS')
    end

end



%% Compile plexon file info, and load strobe event timestamps & values
tic
% General PLX info struct
info = getPlxInfo(plxFilename);

% strobed plexon events loaded as standard component of plxInfo struct
plxStrobeTs     = info.strobes(:,1);
plxStrobeVals   = info.strobes(:,2);

% tic
% % Get strobed plexon events (Omniplex uses 'channel 257' for the strobbed events)
% [~, plxStrobeTs, plxStrobeVals] = plx_event_ts(plxFilename, plx_event_resolve_channel(plxFilename, 'Strobed'));
% 
% if plxStrobeTs == -1
%     fprintf(['\tPlexon MAP system detected. Will use hardcoded channel #257 for strobed sync signals.\n'])
%     [~, plxStrobeTs, plxStrobeVals] = plx_event_ts(plxFilename, 257);
% end
fprintf('\t%2.2f sec to load events from plx file:\t%s\n\n', toc, plxFilename);


%% Reconcile timing among Psychtoolbox(PTB), Plexon(PLX), and Datapixx(DPX)


%% Matchup sets of strobed values to trial "unique numbers" from the PDS & PLX files
% The data in Plexon files often span syncs & data from multiple PLDAPS stimulus files
% Here, we find the unique number values that are sent to Plexon from the PLDAPS computer
% and Spike data that corresponds to this PLDAPS file's trials.
%
% NOTE:  This would all be easier if there was a proper PLDAPS file start & file end sync

% % Unique trial numbers are constructed as follows:      (...from pldapsDefaultTrialFunction.m)
% %         % Construct a unique trial number based on:
% %         %   6 element clock time: [year, month, day, 24hour, minute, second]
% %         unique_number = fix(clock); 
% %         % 	substitute year with trial number (i.e. something actually relevant on the scale of an experimental session)
% %         unique_number(1) = p.trial.pldaps.iTrial; 
% %         %   shift unique numbers into the upper half of our 15-bit strobed word range
% %         unique_number = unique_number + 2^14;
% %         %  ...leaving the lower 16,383 values for easily identifiable event values (p.trial.event).

% Compile PDS trial 'unique_numbers' into a matrix   (sized: ntrial-by-6)
ptbUNum = cell2mat(cellfun(@(X) X.unique_number, PDS.data, 'uni',0)');

% Compile PLX events into the [6] unique number strobe values
uNumStrobes = find(plxStrobeVals>=2^14);    % identify components of unique number strobes by having values >2^14
%         % Hacky fix if need to work OTF:
%         %  uNumStrobes = uNumStrobes(end+1-numel(ptbUNum):end);
if rem(numel(uNumStrobes),6)
    fprintf(2, '/t~!~\nQuantity of unique number strobes is not divisible by 6,\nlikely a crashed file.\n\nAttempting to continue after trimming syncs...\n\t~!~\n');
    jnk = numel(uNumStrobes);
    uNumStrobes = uNumStrobes(1:(jnk-rem(jnk,6)));
end
uNumStrobes = reshape(uNumStrobes, [6,numel(uNumStrobes)/6])';  % reshape into "sixlets" as they were sent
plxUNum = plxStrobeVals(uNumStrobes);

% Convert to datenums corresponding to 8-bit values of each unique_number vector
% (sized: ntrial-by-1)
ptbUDNum = datenum( mod(ptbUNum, 2^8) );
plxUDNum = datenum( mod(plxUNum, 2^8) );

% Match em up (...finally!!!)
[has, plx2ptbUNumIdx] = ismember(ptbUDNum, plxUDNum);
matchedTrialNumbers = find(has);    %plx2ptbUNumIdx(plx2ptbUNumIdx>0);% find(has);


%% %% Compile strobe timestamps from [PTB, DPX, PLX], and compute timestamp conversions
% %  [uNumTimings] becomes matrix unique number timestamps sized n-by-3, as [ptb, datapixx, plx]

uNumTimings = cell2mat(cellfun(@(x) x.datapixx.unique_number_time, PDS.data, 'uni',0)'); % [ptb, datapixx] times (== [nTrials*6, 2])
% Limit PLX syncs to unique numbers from PDS file of interest
%   (i.e. exlude any trial syncs from other files that may be recorded in this same PLX data file)
ii = uNumStrobes(plx2ptbUNumIdx,:)'; 
uNumTimings(:,3) = plxStrobeTs(ii(:));

nt = size(uNumTimings,1);

% Fit conversions between Psychtoolbox(PTB), Datapixx(DPX), & Plexon(PLX) clocks
% using matched unique number timestamps from each device

% Plexon to Psychtoolbox clock
plx2ptbFit = [uNumTimings(:,3), ones(nt,1)] \ uNumTimings(:,1);

% Plexon to Datapixx clock
plx2dpxFit = [uNumTimings(:,3), ones(nt,1)] \ uNumTimings(:,2);


% ----Conversion function handles
% Plexon to PTB
plx2ptb=@(x) x*plx2ptbFit(1) + plx2ptbFit(2);
% PTB to Plexon
ptb2plx=@(x) (x - plx2ptbFit(2))/plx2ptbFit(1);
% Plexon to Datapixx
plx2dpx=@(x) x*plx2dpxFit(1) + plx2dpxFit(2);
% Datapixx to Plexon
dpx2plx=@(x) (x - plx2dpxFit(2))/plx2dpxFit(1);


%% Trim strobes to only those occuring during this file
% i.e.  AFTER PLDAPS experiment start time
%       &
%       BEFORE last datapixxTRIALEND sync**
% 
%   NOTE: ** Adds 1 sec slop time to allow for other modules to complete potential cleanUpAndSave syncs
%   TODO: Would be nice if p.run initiated a true exptEnd sync after all module execution & saving was complete
 
% Convert everything to PTB machine time  (Really?? ...could be utterly simplifying, or utterly damning?)
plxStrobeTs = plx2ptb(plxStrobeTs);

exptStart = PDS.baseParams.session.experimentStart;
exptEnd = PDS.data{end}.timing.datapixxTRIALEND(1) +1; % ** 1 sec slop time added to last trialend timestamp **

keepers = plxStrobeTs>=exptStart & plxStrobeTs<=(exptEnd);
plxStrobeTs(~keepers) = [];
plxStrobeVals(~keepers) = [];


% get TRIALSTART & TRIALEND indices
plxTrialStart = plxStrobeTs(plxStrobeVals==PDS.baseParams.event.TRIALSTART);
ptbTrialStart = cellfun(@(x) x.timing.datapixxTRIALSTART(1), PDS.data)';
dpxTrialStart = cellfun(@(x) x.timing.datapixxTRIALSTART(2), PDS.data)';

plxTrialEnd = plxStrobeTs(plxStrobeVals==PDS.baseParams.event.TRIALEND);


%% plot timing checks if desired
if verbose
    % plot sync drift between devices
    lw = 2;
    figure(999); clf
    hold on
    % Datapixx to Plexon error (in grey;  AWK: must convert Plexon syncs back to native Plexon time for comparison with Datapixx time)
    reconErrorPLX = ( ptb2plx(plxTrialStart) - dpx2plx(dpxTrialStart) )*1000; %ms 
    plot(ptbTrialStart-exptStart, reconErrorPLX, ':', 'color',.3*[1 1 1], 'linewidth',lw+1);
    recErrDP2PL = max(abs(reconErrorPLX));
    L{1} = sprintf('DPX-to-PLX   [%0.4f,  %0.4f] ms',recErrDP2PL, std(reconErrorPLX));
    
    % Plexon to Psychtoolbox/PLDAPS error (in BLUE; this conversion has already been applied)
    reconErrorPTB = ( plxTrialStart - ptbTrialStart )*1000; %ms
    plot(ptbTrialStart-exptStart, reconErrorPTB, 'b-', 'linewidth',lw);
    recErrPTB2PL = max(abs(reconErrorPTB));
    L{2} = sprintf('PLX-to-PTB   [%0.4f,  %0.4f] ms',recErrPTB2PL, std(reconErrorPTB));
    title(PDS.baseParams.session.file, 'interp','none')
    ylabel('Clock sync error (ms)')
    xlabel('Time rel. to Expt Start (sec)')
    box off
    
    lh = legend(L);
    set(lh, 'linewidth',.5, 'location','best')
    try % avoid error on older Matlab graphics engines
        lh.Title.String = 'Recon Error [max, std]';
    end

end

%% Prepare outputs

    % already trimmed to this file duration, just clip out unique numbers
    allSyncIdx = plxStrobeVals<2^14;
    
% .strobe == [ts, value, trial#];     % (...trial# appended below)
nsync.strobe = [plxStrobeTs(allSyncIdx), plxStrobeVals(allSyncIdx)];


% Trial [start, end] timestamps
nsync.plxTrialTs = nan(length(plxTrialStart), 2);
nsync.plxTrialTs(:,1) = plxTrialStart;
nsync.plxTrialTs(1:length(plxTrialEnd), 2) = plxTrialEnd;
% Hacky assignment b/c of occasional missed plxTrial end syncs.
%   --might be due to pauses or strange hiccups, but no reason for it to derrail the whole show

% Create vector of sync trial#
syncTrialNumber = nan(sum(allSyncIdx),1);
% for each trial
for i = 1:size(nsync.plxTrialTs,1)
    % syncs during this trial (this assignment method is redundant, but should be robust to possible missing trialEnd syncs)
    ii = (nsync.strobe(:,1) >= nsync.plxTrialTs(i,1));  % &  (nsync.strobe(:,1) <= nsync.plxTrialTs(i,2));
    syncTrialNumber(ii) = i;
end
% append to strobe matrix
nsync.strobe(:,3) = syncTrialNumber;

%% conversion fxns
nsync.plx2ptb = plx2ptb;
nsync.ptb2plx = ptb2plx;
nsync.plx2dpx = plx2dpx;
nsync.dpx2plx = dpx2plx;
%% syncing info
nsync.info = info;
nsync.info.pdsSrc = pdsFilename;
%% fit data
nsync.info.pl2ptbFit = plx2ptbFit;
nsync.info.pl2dpFit = plx2dpxFit;


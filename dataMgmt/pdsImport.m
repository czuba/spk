function [pds, pdsPath] = pdsImport(filename, theseFields, verbosity)
% function pds = pdsImport(filename, theseFields)
% 
% Make loading .PDS data files (from PLDAPS) more robust & user friendly.
% 
% INPUTS:
%   [filename]  Can be fullpath to PDS file, dir of PDS files to choose from, or standard PLDAPS session directory (/yyyymmdd)
% 
% OUTPUTS:
%   [pds]       Struct of standard PDS data fields
%               - info struct [.info] based on standard PLDAPS content
%               - will attempt to parse additional info based on standard file naming convention:
%                       <subject>_<yyyymmdd><alpha>_<gridLoc>-<depth>_<pdsTime>.PDS
%                 e.g.  tbc_20190131d_M1P3-5200_1520.PDS
%                 ...as would result from executing the following PLDAPS experiment at time 15:20
%                       >> p = modularDemo.doRfPos_gabGrid(["tbc","d_M1P3-5200"]);
% 
%   [pdsPath]   fullfile path
% 
% 2018-10-26  TBC  Wrote it.
% 2019-xx-xx  TBC  Rolling updates...

% TODO: write the gui thing to list available files & show a summary when each file selected
% -- must be version of PDS file that saves as "-struct"
%    to allow indiv. fields to be loaded w/o getting the whole damn thing.
% 
% -- Change to uigetfile.m instead of homerolled chooseFile subfunction...


%% Parse inputs
if nargin <1
    filename = [];
end

if nargin <2
    theseFields = [];
end
if nargin <3 || isempty(verbosity)
    verbosity = 1;
end

if isempty(filename)
    pdsPath = pwd;
    % Each recording day has its own "pds" directory auto-generated
    if exist( fullfile(pdsPath, 'pds'), 'dir')
        pdsPath = fullfile(pdsPath, 'pds');
    end
    
elseif isfolder(filename)
    % directory specified, select w/gui
    pdsPath = filename;
    filename = [];
    
elseif isfile(fullfile(pwd, 'pds', filename))
    % standard pds location w/in recording session dir
    pdsPath = fullfile(pwd, 'pds');
    
elseif isfile(filename)
    % fully specified file, or already in pds directory
    if contains(filename, filesep) % is a full or partial path
        [pdsPath, filename, ext] = fileparts(filename);
        filename = [filename, ext];
    else
        pdsPath = pwd;
    end
    
else 
    % Don't know...
    keyboard
end


if isempty(filename)
    fd = dir( fullfile(pdsPath, '*.PDS') );
else
    fd = dir( fullfile(pdsPath, filename) );
    
end

% Load PDS file
if ~isempty(fd)
    if length(fd)>1
        [fd, ok] = chooseFile(fd, 'Select PDS source file:');
        if ~ok
            fprintf(2, 'No PDS files found in %s, or user canceled selection\n', pdsPath);
            return
        end
    end
    
    switch lower(fd.name(end-2:end))
        case 'mat'
            % use pre-compiled data struct
            % ...this method not really fully coded, but might "just work"
            thisfile = fullfile(pdsPath, [fd.name(1:end-3),'mat']);
            if verbosity
                fprintf('Loading pre-compiled analysis struct:\n\t%s\n', thisfile);
            end
            pds = load( thisfile );
            
        case 'pds'
            % load PDS struct
            pdsPath = fullfile(pdsPath, fd.name);
            if verbosity
                fprintf('Loading PDS file:\n\t%s\n', pdsPath);
            end
            
            % Must be modern PDS file (>=glDraw branch). Damn the torpedos!!
            if isempty(theseFields)
                pds = load(pdsPath, '-mat');
                if length(fieldnames(pds))==1
                    % likely temp file or not saved correctly...try unpacking
                    fn = fieldnames(pds);
                    pds = pds.(fn{1});
                end
            else
                pds = load(pdsPath, '-mat', theseFields);
            end
            if isa(pds,'pldaps')
                return
                % file saved without converting to struct output (...only expected for "trial00000" TEMP file)
                pds = pds.save;
                
            end
            
            try
                % PDS stimulus file info
                if isfield(pds,'info')
                    info = pds.info;
                else
                    info = struct;
                end
                [~, info.pdsName] = fileparts(pds.baseParams.session.file);
                
                info.stimType = {pds.baseParams.pldaps.modNames.currentStim{1}, pds.baseParams.session.caller.name};
                info.viewDist = pds.baseParams.display.viewdist;
                
                % default stereo probe depths
                trodeDepth = fliplr(kron(0:100:1500, [1,1]))';
                
                % try to expand file name components (may fail)
                fnParts = strsplit(info.pdsName, {'_','-'});
                for i = 1:length(fnParts)
                    switch i
                        case 1
                            info.subj = fnParts{i};
                        case 2
                            info.sesh = fnParts{i};
                        case 3
                            if contains(fnParts{i}, {'M','P','L','A'})
                                info.gridLoc = fnParts{i};
                                isOffset = -0.5*contains(info.gridLoc,'o');
                                info.gridXY = nan([1,2]);
                                % ML(X) grid loc in mm
                                if contains(info.gridLoc, 'M')
                                    ii = strfind(info.gridLoc,'M')+1;
                                    info.gridXY(1) = isOffset + str2num(info.gridLoc(ii));
                                elseif contains(info.gridLoc, 'L')
                                    ii = strfind(info.gridLoc,'L')+1;
                                    info.gridXY(1) = -(isOffset + str2num(info.gridLoc(ii)));
                                end
                                
                                % AP(Y) grid loc in mm
                                if contains(info.gridLoc, 'A')
                                    ii = strfind(info.gridLoc,'A')+1;
                                    info.gridXY(2) = isOffset + str2num(info.gridLoc(ii));
                                elseif contains(info.gridLoc, 'P')
                                    ii = strfind(info.gridLoc,'P')+1;
                                    info.gridXY(2) = -(isOffset + str2num(info.gridLoc(ii)));
                                end
                                
                            end
                            
                        case length(fnParts)
                            % last is always PDS time
                            info.pdsTime = fnParts{i};
                        case 4
                            % assume depth of probe
                            % strip "hm" handmapping designation
                            dstr = fnParts{i};
                            ii = [strfind(dstr,'h'), strfind(dstr,'m')];
                            dstr(ii) = [];
                            try
                                info.siteDepth = str2num(dstr);
                                info.trodeDepth = info.siteDepth-trodeDepth;
                            end
                    end
                end
            end
            pds.info = info;
            
            if verbosity
                fprintf('\tDone.\n')
            end
    end
    
end


end %main function

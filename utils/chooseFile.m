function [fd, ok] = chooseFile(fd, titl, selMode, defSel)
% function [fd, ok] = chooseFile(fd, titl, selMode, defSel)
% 
% Simple file selection dialog to select outputs returned by dir()
%   fd == file struct output from:  fd = dir(____);
%   titl == file selection window title/prompt [def: 'Choose file (selMode):']
%   selMode = 'single' or 'multiple'[default]
%   defSel = default selection index(es) or 'all'; [last]
%
% Outputs:
% [fd] dir input struct trimmed to selection. (empty if canceled)
% [ok] 1 if pressed "ok" or return key, 0 if pressed "cancel" or esc key
% 
% 2018-12-00  TBC  Functionified it from many prev. subfunction uses

if nargin<3 || isempty(selMode) || strcmpi(selMode,'multiple')
    selMode = 'multiple';
else
    selMode = 'single';
end

if nargin<2 || isempty(titl)
    titl = sprintf('Choose file (%s):',selMode);
end

if nargin<4 || isempty(defSel)
    % select the last/most recent one by default
    defSel = length(fd);
elseif strcmp(defSel,'all')
    defSel = 1:length(fd);
elseif islogical(defSel)
    defSel = find(defSel);
end


winsz = [400, 60+12*length(fd)];
winsz(winsz>400) = 600;
[fdsel, ok] = listdlg('liststring',{fd.name}, 'selectionmode',selMode, 'initialValue',defSel, 'listsize',winsz, 'Name', titl);
if ok && ~isempty(fdsel)
    fd = fd(fdsel);
else
    disp('User canceled selection')
    %         % use most recent file match/sorting
    %         fd = fd(end);
end


end

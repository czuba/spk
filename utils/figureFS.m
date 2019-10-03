function hout = figureFS(h, ori, type)

% function h = figureFS(h, ori, type)
%   h       = handle to figure (non destructive). Opens new new figure by default.
%   ori     = paper orientation:  'portrait' or 'landscape'(def)
%   type    = paper type: string 'A4', 'B5', 'usletter', 'tabloid'(def), etc.
%               (...or width x height in inches)
% 
% Adjust figure settings so that print/save versions will default to fill the page.
% 
% 2012-12-09  TBC  Wrote it. (Finally!)
% 2014-07-16  TBC  Updates & input info


if ~exist('h','var')
    h = figure;
elseif isempty(h)
    h = figure;
elseif iscell(h)
    if isempty(h{1})
        h{1} = figure;
    else
        figure(h{1});
    end
    set(h{1}, 'name',h{2});
    h = h{1};
elseif ischar(h)
    nm = h;
    h = figure;
    set(h, 'name',nm);
else
    figure(h)
end

if ~exist('ori','var') || isempty(ori)
    ori = 'landscape';
end

if ~exist('type','var') || isempty(type)
    type = 'tabloid';
end


set(h,'PaperOrientation', ori);
    orient(ori); %...and again (thanks Matlab)
    
try
    set(h,'PaperType', type);  %'usletter');%
catch
    set(h, 'PaperType','<custom>', 'PaperSize', type)
end


set(h,'PaperUnits', 'inches')
sz = get(h,'PaperSize');
set(h,'PaperPositionMode','manual');
n = min([.1*sz, 0.5]);
set(h,'PaperPosition',[n/2 n/2 sz(1)-n sz(2)-n]);

% just keep doing this till it sticks
set(h,'PaperOrientation', ori);
orient(ori); %...and again (thanks Matlab)

% set(h,'PaperUnits', 'normalized')


% silence unneeded outputs
if nargout>0
    hout = h;
end

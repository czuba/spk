function [spx, spy, figsz, spMargin, pbaxr, spbase, figHeight, figWidth] = mkProbeSubplots(nsp, oriFlag)
% function [spx, spy, figsz, spMargin, pbaxr, spbase, figHeight, figWidth] = mkProbeSubplots(nsp, oriFlag)
% 
% A scrappy attempt to standardize probe data subplotting
% - initializes set of subplot size & spacing parameters for use with subplot_tight.m
% 
% INPUTS:
%   [nsp]       Total number of subplots desired (usually nsp==nunits)
%   [oriFlag]   Optional second input flag:
%       'g' (default) produces 8-by-n grid of subplots
%       'v' for 2-by-n vertical pairs; arranged like unsorted channels of stereoprobe
%       'h' for n-by-2 rows; marginally more 'monitor-friendly', primarily used back in the <=16 channel days)
% 
% OUTPUTS:
%   [spx, spy]  Number of subplots in x & y directions
%               - *NOTE* quirky order of subplot.m inputs means these should be used as:  subplot(spy, spx,...)
%   [figsz]     Size of overall figure in inches, as [widht, height]
%               - use with figureFS.m helper function(...figure 'F'or 'S'aving)
%   [spMargin]  total normalized figure space apportioned to gaps between subplots (normalized units; def=[0.03 0.025])
%               - passed as last/4th input to subplot_tight.m
%   [pbaxr]     Plot box aspect ratio, fixed default=[4,3,1] is appropriate for most tuning plots
%   ----outdated/legacy----
%   [spbase]    base subplot # offset, used for 2-row 'h' plotting with related axes above one another (e.g. RFposition & direction tuning)
%   [figHeight, figWidth]   ...just components of [figsz]
% 
% 
% EXAMPLE USAGE:
% %----
%   [spx, spy, figsz, spMargin, pbaxr, spbase] = mkProbeSubplots(nunits);
%   H1 = figureFS(99, 'portrait', figsz);
%   ...
%   for u = 1:nunits
%       figure(H1);
%       sp = subplot_tight(spy, spx, u+spbase, spMargin);
% 
%       ... % do your fitting & plotting
% 
%       set(sp,'plotboxaspectratio', pbaxr);
%   end
%   saveFigTriplet(1, 'Extra info/description/source/etc', {'png','pdf});
% %----
% 
%   
% see also: figureFS, subplot_tight, addWfInset, plotUnitDirTune
% 
% 2021-xx-xx  TBC  Wrote it.
% 


if nargin<2 || isempty('oriFlag')
    oriFlag = 'g'; % standard 8-by-n grid of subplots
elseif ischar(oriFlag)
    oriFlag = lower(oriFlag);
else
    errstr = sprintf('oriFlag must be one of the following:\n\t''h'' (n-by-2 horizontal)\n\t''v'' (2-by-n vertical)\n\t''g'' (8-by-n grid)');
    error(errstr);
end

% subplot arrangements
spx = nsp;
% if isempty(findobj(H, 'type','axes','tag',sprintf('ch%drf',1)))
%     spbase = 0;
%     spx = spx/2;
% else
%     spx = nunits;
    spbase = 0;
% end


if strcmpi(oriFlag,'g') || (spbase==0 && nsp~=32)
    spx = 8;
    spy = ceil(nsp/spx);
elseif strcmpi(oriFlag,'v')
    spx = 2;
    spy = ceil(nsp/spx);
else %if strcmpi(oriFlag,'h')
    spy = 2;
    spx = ceil(nsp/spy);
end

% plotbox aspect ratio
% pbaxr = [3,2,1];
pbaxr = [4,3,1];

figHeight = pbaxr(1)*spy;
figWidth = pbaxr(2)*spx + 2;
figsz = [figWidth, figHeight];

spMargin = [.03 .025];%1/nunits


end %main function

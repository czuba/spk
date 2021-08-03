function [spx, spy, figsz, spMargin, pbaxr, spbase, figHeight, figWidth] = mkProbeSubplots(nunits, oriFlag)
% function [spx, spy, figsz, spMargin, pbaxr, spbase, figHeight, figWidth] = mkProbeSubplots(nunits, oriFlag)
% 
% scrappy attempt to standardize probe data subplotting
%
% USAGE:
% % ...for u = 1:nunits
% % ...
%      sp = subplot_tight(spy, spx, u+spbase, spMargin);
%
%

if nargin<2 || isempty('oriFlag')
    oriFlag = 'g';
elseif ischar(oriFlag)
    oriFlag = lower(oriFlag);
else
    errstr = sprintf('oriFlag must be one of the following:\n\t''h'' (n-by-2 horizontal)\n\t''v'' (2-by-n vertical)\n\t''g'' (8-by-n grid)');
    error(errstr);
end

% subplot arrangements
spx = nunits;
% if isempty(findobj(H, 'type','axes','tag',sprintf('ch%drf',1)))
%     spbase = 0;
%     spx = spx/2;
% else
%     spx = nunits;
    spbase = 0;
% end


if strcmpi(oriFlag,'g') || (spbase==0 && nunits~=32)
    spx = 8;
    spy = ceil(nunits/spx);
elseif strcmpi(oriFlag,'v')
    spx = 2;
    spy = ceil(nunits/spx);
else %if strcmpi(oriFlag,'h')
    spy = 2;
    spx = ceil(nunits/spy);
end

% plotbox aspect ratio
% pbaxr = [3,2,1];
pbaxr = [4,3,1];

figHeight = pbaxr(1)*spy;
figWidth = pbaxr(2)*spx + 2;
figsz = [figWidth, figHeight];

spMargin = [.03 .025];%1/nunits

function hout = xyunity(h, xs, ys)

if nargin<1 || isempty(h) || ~ishandle(h)
    h = gca;
end

hld = get(h, 'NextPlot');
hold on

col = .75*[1 1 1];
len = 1e5;

x = len*[-1,1;0,0;-1,1;eps,1]';
y = len*[0,0;-1,1;-1,1;eps,1]';

if nargin==3
    if isscalar(xs) && isscalar(ys)
        % assume scalar shift of origin
        x = x+xs;
        y = y+ys;
    else
        x = [x,xs];
        y = [y,ys];
    end

elseif nargin==2
    % assume xs is a scalar shift origin
    x = x+xs;
    y = y+xs;
end

hl = plot(h, x,   y, '--', 'linewidth',.5, 'color',col, 'xliminclude','off', 'yliminclude','off');

% if nargin<=1
%     hl = plot(h, len*[-1,1;0,0;-1,1;eps,1]',   len*[0,0;-1,1;-1,1;eps,1]', '--', 'linewidth',.5, 'color',col, 'xliminclude','off', 'yliminclude','off');
% else
%     hl = plot(h, [len*[-1,1;0,0;-1,1;eps,1]', xs],   [len*[0,0;-1,1;-1,1;eps,1]', ys], '--', 'linewidth',.5, 'color',col, 'xliminclude','off', 'yliminclude','off');
% end

% return hold state as before
set(h, 'NextPlot',hld);

if nargout
    hout = hl;
end

return

function fmtPlot(H, type)
% function fmtPlot(H, type)
%
% Tighten up & standardize axis formatting
% [type]
% 
% 2020-09-24 TBC Wrote it.
% 

if nargin<1 || isempty(H)
    H = gca;
elseif nargin<2 && ischar(H)
    type = H;
    H = gca;
end

if ~exist('type','var') || isempty(type)
    type = 'small';
end


switch lower(type)
    case {'small','xsmall'}
        %if strcmp(H.XLimMode,'auto') && strcmp(H.YLimMode,'auto')
        axis(H,'tight')
        %end        
        fsz = 10;
        axgap = 0.05;
        if strcmpi(type(1),'x')
            % xtra small
            fsz = 8;
            axgap = 0.01;
        end
        set(H, 'fontsize',fsz, 'linewidth',.65, 'box','off')
        % Set xy limits
        yl = ylim;
        xl = xlim;
        if ~isempty(findobj(gca,'type','image'))
            % skip if image
        else
            ylim(yl + axgap*diff(yl)*[-1,1]);
            xlim(xl + axgap*diff(xl)*[-1,1]);
        end
end

end %main function

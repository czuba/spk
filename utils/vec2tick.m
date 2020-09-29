function st = vec2tick(vec, fmt)
% function st = vec2tick(vec, fmt)
% 
% convert vector[vec] of numbers to a cell of strings.
%   [fmt] = text format string (def= '%2.2f ')
% only outputs single row of stringified numbers (no matrix possible...)
% 
% 2017-01-04  TBC  Wrote it.

if nargin<2 || isempty('fmt') || ~ischar(fmt)
    % must be a formatting string
    fmt = '%.5g ';
elseif (fmt(end)*1) ~= (1*' ')
    % format must end in space
    fmt = [fmt, ' '];
end
    
st = textscan( num2str(vec(:)', fmt), '%s');
st = st{:};


% % % st = strread( num2str(vec, fmt), '%s');  % orig, but matlab recommended against strread.m (...dont know why)

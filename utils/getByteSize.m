function b = getByteSize(thisVar, returnUnits, returnStr)         %#ok<INUSL>
% function b = getByteSize(theVariable, returnUnits, fid)
%
% getByteSize returns the mem.usage of input [thisVar]
%
% INPUTS
%   [thisVar]     input variable
%   [returnUnits] string input determining byte units to return
%       - ('byte', 'kb', 'mb', 'gb', 'tb')
%       - if empty, will autoselect appropriate scaled units *for each *return value (not all same)
%   [returnStr]   logical flag to return a string (incl units) instead of numerical values (def=0)
%       - string output automatically pads to 8 
% 
% Output is written to screen if no output arguments, else values [only] returned in [b]
% - ...if no untis provided, [b] return argument units will be ambiguous
%
% 2021-07-xx  TBC  copied & revised from Matlab File Exchange inline solution
% 

if nargin<3 || isempty(returnStr)
    returnStr = 0;
end

s = whos('thisVar');
b = s.bytes;
if nargin<2 || isempty(returnUnits)
    scale = floor(log(b)/log(1024));
    switch scale
        case 0
            returnUnits = ' b';
        case 1
            returnUnits = 'kb';
        case 2
            returnUnits = 'mb';
        case 3
            returnUnits = 'gb';
        case 4
            returnUnits = 'tb';
        case -inf
            % Size occasionally returned as zero (eg some Java objects).
            returnUnits = ' b';
            % warning('Size occasionally returned as zero (eg some Java objects). Bytes assumed');
        otherwise
            returnUnits = 'pb'; %petabytes
            warning('Over 1024 petabyte. petabytes assumed');
    end
end

switch returnUnits
    case {' b','byte','bytes'}
        b = s.bytes;
    case {'kb','kbs','kilobyte','kilobytes'}
        b = b / 1024;
    case {'mb','mbs','megabyte','megabytes'}
        b = b / 1024^2;
    case {'gb','gbs','gigabyte','gigabytes'}
        b = b / 1024^3;
    case {'tb','tbs','terabyte','terabytes'}
        b = b / 1024^4;
    case {'pb','pbs','petabyte','petabytes'}
        b = b / 1024^5;
    otherwise
        returnUnits = 'bytes';
end
if ~nargout
    % ONLY print to screen if no outputs
    fprintf('%s %s\n', num2str(b), returnUnits);
    clear b
    return
elseif returnStr
    % transform to string output
    b = sprintf('%s %s', pad(num2str(b, '%8.3g'),8), returnUnits);
end

end % main function

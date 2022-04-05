function [dirbias,dirpref,oribias,oripref,quadbias,quadpref] = orivecfit(ori,response,baseline)
%
% function [dirbias,dirpref,oribias,oripref,quadbias,quadpref] = 
%   orivecfit(ori,response,baseline)
% OR
% function [dirbias,dirpref,oribias,oripref,quadbias,quadpref] = 
%   orivecfit(ori,response)
%
% dirpref, oripref, and quadpref are vector fits of direction
%  on a 360, 180, or 90 degree scale, respectively
% dirbias, oribias, and quadbias are measures of the selectivity
%  of those vector fits
%
% ori and response are 1D vectors containing the orientations in
% degrees and the response firing rates. If absent, the baseline
% vector is assumed to be zero (useful for F1 data). Otherwise, it
% can be a vector of the baseline responses or a single number.
%

% Matthew A. Smith
% Revised: 20011101

if (nargin < 3) % F1 data have no baselines
  baseline = 0;
end

if ((length(baseline)) > 1)
  baseline = mean(baseline);
end

oriR  = ori .* (pi/180);

RDnum = 0;
IDnum = 0;
ROnum = 0;
IOnum = 0;
RQnum = 0;
IQnum = 0;
denom = 0;

for I = 1:length(response)
  y = cos(oriR(I));
  x = sin(oriR(I));
  w = cos(2*oriR(I));
  v = sin(2*oriR(I));
  p = cos(4*oriR(I));
  q = sin(4*oriR(I));
  RDnum = RDnum + (response(I) - baseline)* y;
  IDnum = IDnum + (response(I) - baseline)* x;
  ROnum = ROnum + (response(I) - baseline)* w;
  IOnum = IOnum + (response(I) - baseline)* v;
  RQnum = RQnum + (response(I) - baseline)* p;
  IQnum = IQnum + (response(I) - baseline)* q;
  denom = denom + abs(response(I) - baseline);
end

Realdir = RDnum/denom;
Imdir = IDnum/denom;
dirbias = sqrt(Realdir^2 + Imdir^2);
dirpref = (180/pi * atan2(Imdir, Realdir));
if (dirpref < 0)
  dirpref = dirpref + 360;
end

% ...be speedy
if nargout<3
    return
else
    Realori = ROnum/denom;
    Imori = IOnum/denom;
    oribias = sqrt(Realori^2 + Imori^2);
    oripref = .5 * (180/pi * atan2(Imori, Realori));
%     if (oripref < 0)
%         oripref = oripref + 360;
%     end
    % ...still more speedy
    if nargout<5
        return
    else
        
        Realq = RQnum/denom;
        Imq = IQnum/denom;
        quadbias = sqrt(Realq^2 + Imq^2);
        quadpref = .25 * (180/pi * atan2(Imq, Realq));
        if (quadpref < 0)
            quadpref = quadpref + 360;
        end
    end
end
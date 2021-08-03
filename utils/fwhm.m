function width = fwhm(x,y)

% function width = fwhm(x,y)
%
% Full-Width at Half-Maximum (FWHM) of the waveform y(x)
% and its polarity.
% The FWHM result in 'width' will be in units of 'x'
%
%
% Rev 1.2, April 2006 (Patrick Egan)
% 2019-05  TBC  Partially clean up this [really!] scrappy fxn grepped from matlab file exchange

%%
width = NaN;

y = normalize(y);   %y / max(y);
N = length(y);
lev50 = 0.5;

[~, centerindex] = maxk(y,1);
i = 2;
while sign(y(i)-lev50) == sign(y(i-1)-lev50)
    i = i+1;
end
%first crossing is between v(i-1) & v(i)
halfstp = (lev50-y(i-1)) / (y(i)-y(i-1));
tlead = x(i-1) + halfstp*(x(i)-x(i-1));

%start search for next crossing at center
i = centerindex+1;                    
while ((sign(y(i)-lev50) == sign(y(i-1)-lev50)) & (i <= N-1))
    i = i+1;
end

if i == N
    fprintf(2, 'Step-Like Pulse, no second edge')
    width = nan;
    return
end

%     disp('Pulse is Impulse or Rectangular with 2 edges')
halfstp = (lev50-y(i-1)) / (y(i)-y(i-1));
ttrail = x(i-1) + halfstp*(x(i)-x(i-1));
width = ttrail - tlead;


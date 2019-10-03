function f = vmDouble(x,xdata)
% x     = von mises parameters (can be 5-by-1 or 5-by-n if passing multiple parameter sets)
% xdata = xvalues to evaluate function at (radians)
%
%       params are: [m, a1, b, t, a2]
%             [mean, majorAmp, spread, thetaPref(rad), minorAmp]
%
if min(size(x))==1   % careful vectorizing
    m = x(1);
    a1 = x(2);
    b = x(3);
    ThetaP = x(4);
    a2 = x(5);
    f = m + a1*exp(b*(cos(xdata-ThetaP)-1))+a2*exp(b*(cos(xdata-ThetaP+pi)-1));
    
else
    
    m = x(1,:);
    a1 = x(2,:);
    b = x(3,:);
    ThetaP = x(4,:);
    a2 = x(5,:);
    
    f = nan(length(m), length(xdata));
    for i = 1:length(m)
        f(i,:) = m(i) + a1(i)*exp(b(i)*(cos(xdata-ThetaP(i))-1))+a2(i)*exp(b(i)*(cos(xdata-ThetaP(i)+pi)-1));
    end
end
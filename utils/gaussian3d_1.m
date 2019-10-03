function f = gaussian3d_1(x,xdata)
% 3D gaussian fxn [x,y,response] with COMMON sigma for x & y dimensions
%

th = x(1);  % baseline
a =  x(2);  % amplitude
x0 = x(3);  % x-mean
y0 = x(4);  % y-mean
sd = x(5);  % sigma
%sdy = x(6);
%th = 0;

f = th + a*exp(-.5*(xdata{1}-x0).^2/sd - .5*(xdata{2}-y0).^2/sd);

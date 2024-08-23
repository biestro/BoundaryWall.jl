function [x,y,xm,ym] = createParabolicBilliard(xi0,eta0,N)
% curve 1
xi = xi0; eta = eta0 * (2/N) * (0: N/2  ); % div entre 4 para que sean N/2 puntos
x = (xi.^2 - eta.^2)/2; y = xi.*eta;


% curve 2
xi = xi0 * (2/N) * (0:N/2); eta = eta0;
x1 = (xi.^2 - eta.^2)/2; y1 = xi.*eta;
% append parabolas
x = [x x1(end-1:-1:1)]; y = [y y1(end-1:-1:1)];
x = [x flip(x)];        y = [y -y];
close all
% interpolate for equal arc lengths
[x, y]  = subdivide(x,y, 'spline');
[xm,ym] = midpoints(x,y);
%[xi_m, eta_m] = cart2par(xm,ym);
end

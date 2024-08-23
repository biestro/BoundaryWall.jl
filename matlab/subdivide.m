
function [xs, ys] = subdivide(X,Y, method)
dxdt = gradient(X, X(2) - X(1));
dydt = gradient(Y, X(2) - X(1));
s = cumtrapz( sqrt(dxdt.^2 + dydt.^2 ) );

S = transpose(s);
V = [ transpose(X), transpose(Y) ];

N = length(s);
Xq = linspace(s(1),s(end),N); % equally spaced indices

Vq = interp1(S,V,Xq, method); 

xs = Vq(:,1);
ys = Vq(:,2);
end
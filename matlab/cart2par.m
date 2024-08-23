function [xi, eta] = cart2par(x,y)
% convert cartesian to parabolic coordinates
% convention to have the negative x axis as the 
z = x + 1j * y; % convert to imaginary unit
z = sqrt(2*z);
xi = real(z);
eta = imag(z);
end
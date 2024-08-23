%%   DOMAIN   


%%   CONFOCAL BOUNDARY
xi0 = 3.0; eta0 = 2.0;
[x,y,xm,ym] = createParabolicBilliard(xi0, eta0, 400);

%% EIGENVALUES (Villarreal, et al.)

eigvals =[
    0.403, 0.977;
    0.805, 1.877;
    1.194 ,-0.116;
    1.365, 2.870;
    1.798, 0.107;
    2.064, 3.548;
    2.405 ,-0.531;
    2.629, 1.003;
    2.900, 4.067;
    3.181 ,-0.597;
    3.586, 1.557;
    3.872, 4.475;
    4.004 ,-0.867;
    4.320, 0.316;
    4.697, 2.066;
    4.947 ,-1.004;
    4.976, 4.807;
    5.491, 0.640;
    5.953, 2.498;
    6.213, 5.082;
    6.463, 0.009;
    6.857, 1.076;
    7.350, 2.869;
    7.806, 0.130;
    8.368, 1.452;
    8.887, 3.191;
    9.432, 0.525;
    10.03, 1.796;
    11.17, 0.821;
    11.83, 2.106;
    13.08, 1.123;
    15.15, 1.402;];

k = sqrt(eigvals(:,1)); % obtain k
a = eigvals(:,2) .* k/2; % obtain a
eig_ka = [k, a]; % useful numbers, wavenumber and constant a

%%   SCARRING

% constants
M_scarring    = 2;
PHI0_scarring = pi/2;
p_winding     = 2;
q_turning     = 5;


m0 = 10; n0 = eigvals(5,:);
k_m0n0 = besselzero(m0,n0,1);

% domain
nx = 450; ny = 400;
ydom = linspace(-6.5,6.5,nx);
xdom = linspace(-3,5,ny);

[X,Y] = meshgrid(xdom,ydom);
X = X(:);  Y = Y(:);
[XI,ETA] = cart2par(X,Y);

PSI = zeros(size(X));

eigenfun = @(k_,a_,xi_,eta_) Ue(xi_,eta_,sqrt(2*k_),a_) + 1j*Uo(xi_,eta_,sqrt(2*k_),a_);


i_ = 1;
for K = -M_scarring:M_scarring
    m0_prime = m0 + K*q_turning;
    n0_prime = n0 - K*p_winding;
    %k_mn = besselzero(m0_prime, n0_prime,1)/circle_radius;
    %k_mn = k_mn(end);
    k_mn = eig_ka(m0_prime, 1);
    a_mn = eig_ka(m0_prime, 2);
    PSI = PSI + exp(1j*K*q_turning*PHI0_scarring) .* eigenfun(m0_prime, k_mn, R,TH) / (sqrt(pi)* circle_radius * besselj(m0_prime-1,k_mn*circle_radius));
    % PSI(:,:,i_) = exp(1j*K*q_turning*PHI0_scarring) * 
    i_ = i_ + 1;
end



inside = R < circle_radius;
close all
surf(xdom, ydom, abs(reshape(PSI .* inside, nx, ny)))
colormap(crameri("grayC"))
hold on
plot3(x,y,ones(length(x),1),"k")
view(2);
shading interp;
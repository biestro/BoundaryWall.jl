%%   DOMAIN
ydom = linspace(-6.5,6.5,130);
xdom = linspace(-3,5,130);

[X,Y] = meshgrid(xdom,ydom);
[XI,ETA] = cart2par(X,Y);

%%   CONFOCAL BOUNDARY
xi0 = 3.0; eta0 = 2.0;
[x,y,xm,ym] = createParabolicBilliard(xi0, eta0, 400);



%%   EIGENSTATE

beta    = 5.082;
wavenum = sqrt(6.213);
a       = beta*wavenum/2;
SIGMA   = sqrt(2*wavenum);

wave    = Ue(XI,ETA,SIGMA,a);

outside  = (XI < xi0) .* (abs(ETA) < eta0);

%%   PLOTS
clc; close all;
surf(xdom,ydom, abs2(wave) .* outside);
view(2);
hold on
plot3(xm,ym,ones(length(xm),1), "k")
colormap(crameri("grayC"))
axis equal tight;
shading interp;
hold off

%%   SCARRING

% constants
circle_radius = 1.0;
M_scarring    = 2;
PHI0_scarring = pi/2;
p_winding     = 1;
q_turning     = 4;


%m0 = 300; n0 = 100;
m0 = 100; n0 = 5;
k_m0n0 = besselzero(m0,n0,1);

% domain
nx = 450; ny = 400;
ydom = linspace(-circle_radius,circle_radius,nx);
xdom = linspace(-circle_radius,circle_radius,ny);

[X,Y] = meshgrid(xdom,ydom);

X = X(:);
Y = Y(:);

[TH,R] = cart2pol(X,Y);



theta = linspace(0,2*pi,100);
x = circle_radius*cos(theta);
y = circle_radius*sin(theta);


PSI = zeros(size(X));

eigenfun = @(m,k,r,t) besselj(m,k*r) .* exp(1j*m*t);


i_ = 1;
for K = -M_scarring:M_scarring
    m0_prime = m0 + K*q_turning;
    n0_prime = n0 - K*p_winding;
    k_mn = besselzero(m0_prime, n0_prime,1)/circle_radius;
    k_mn = k_mn(end);
    PSI = PSI + exp(1j*K*q_turning*PHI0_scarring) .* eigenfun(m0_prime, k_mn, R,TH) / (sqrt(pi)* circle_radius * besselj(m0_prime-1,k_mn*circle_radius));
    % PSI(:,:,i_) = exp(1j*K*q_turning*PHI0_scarring) *
    i_ = i_ + 1;
end



outside= R >= circle_radius;
PSI(outside) = 0.0;
close all
surf(xdom, ydom, abs2(reshape(PSI, nx, ny)))
colormap(crameri("grayC"))
hold on
plot3(x,y,ones(length(x),1),"k")
view(2);
shading interp;
axis equal tight;

%%
[gradFx, gradFy] = gradient(reshape(PSI, nx, ny), xdom(2)-xdom(1), ydom(2)-ydom(1));
flowX = imag(reshape(PSI, nx, ny) .* gradFx);
flowY = imag(reshape(PSI, nx, ny) .* gradFy);


d = 10;
close all;
quiver(xdom(1:d:end), ydom(1:d:end), flowX(1:d:end,1:d:end), flowY(1:d:end,1:d:end))
hold on
plot(x,y,"k")
axis equal tight;
hold off





%% EXPORT DATA


% create
%h5create(string(datetime("now"))+'_real.h5', '/data', size(PSI));
%h5create(string(datetime("now"))+'_imag.h5', '/data', size(PSI));



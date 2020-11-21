%% APP
close all
clear all
clc

% Valeurs à D0 = 50
d0_50 = 50;
var_50 = 4;
phi0_50 = 15; % degrés

%Valeurs à D0 = 100
d0_100 = 100;
var_100 = 16;
phi0_100 = 30; % degrés

% 1)-------------------------------------------------
% a)------------------------
N = 10000;
theta = rand(1, N) * 2 * pi;

% b)------------------------
figure(1)
histogram(theta)

% c)------------------------
r = (0:0.01:16);
index = 0;
figure
for variance = [0.25 1 4 9 16]
   index = index + 1;
   fR = r / variance .* exp(-r.^2/(2*variance));
   subplot(5, 1, index)
   plot(r, fR)
   subplot_title = strcat('Rayleigh variance = ', num2str(variance));
   title(subplot_title)
end

% d)------------------------
% CDF
cdf_50 = 1- exp(-1 * r.^2 / (2*var_50));
cdf_100 = 1- exp(-1 * r.^2 / (2*var_100));

figure
subplot(2,1,1)
plot(cdf_50)
title('CDF D0 = 50')
subplot(2,1,2)
plot(cdf_100)
title('CDF D0 = 100')

% e)------------------------
dN = 1/N;
p = (0:dN:1-dN);

% CDF Inverse
r_50 = sqrt(-2*var_50*log(1-p));
r_100 = sqrt(-2*var_100*log(1-p));

figure
subplot(2,1,1)
plot(p, r_50)
title('CDF Inverse D0 = 50')
subplot(2,1,2)
plot(p, r_100)
title('CDF Inverse D0 = 100')

% Générateur aléatoire
distribution = rand(1,N);

cloche_50 = [];
cloche_100 = [];
index = 1;
for val = distribution
   val = floor(val * N)+1;
   cloche_50(index) = r_50(val);
   cloche_100(index) = r_100(val);
   index = index+1;
end

% f)------------------------
figure
subplot(2,1,1)
histogram(cloche_50)
title('Histogramme aléatoire D0 = 50')
xlim([-1 18])
subplot(2,1,2)
histogram(cloche_100)
title('Histogramme aléatoire D0 = 100')
xlim([-1 18])

% g)------------------------
figure
subplot(2,1,1)
Rayleigh_th_50 = raylrnd(mean(r_50)*sqrt(2/pi), N, 1);
histogram(Rayleigh_th_50)
title('Rayleigh théorique D0 = 50')
xlim([-1 18])
subplot(2,1,2)
Rayleigh_th_100 = raylrnd(mean(r_100)*sqrt(2/pi), N, 1);
histogram(Rayleigh_th_100)
title('Rayleigh théorique D0 = 100')
xlim([-1 18])
%----------------------------------------------------

%2) -------------------------------------------------
figure
subplot(2,1,1)
scatter(r_50, theta, 1)
title('Couple r et theta avec variance = 4')
subplot(2,1,2)
scatter(r_100, theta, 1)
title('Couple r et theta avec variance = 16')
%----------------------------------------------------

%3) -------------------------------------------------
% % D0 = 50
% d_50 = sqrt(d0_50^2 + 2*d0_50*r_50.*cos(theta) + r_50.^2);
% phi_50 = (phi0_50 + atan(r_50.*sin(theta)./(d0_50 + r_50.*cos(theta)))) / 180 * pi;
% 
% % D0 = 100
% d_100 = sqrt(d0_100^2 + 2*d0_100*r_100.*cos(theta) + r_100.^2);
% phi_100 = (phi0_100 + atan(r_100.*sin(theta)./(d0_100 + r_100.*cos(theta)))) / 180 * pi;

% D0 = 50
d_50 = d0_50 + r_50.*cos(theta);
phi_50 = (phi0_50 + r_50.*sin(theta))/ 180 * pi;

% D0 = 100
d_100 = d0_100 + r_100.*cos(theta);
phi_100 = (phi0_100 + r_100.*sin(theta))/ 180 * pi;

%----------------------------------------------------

%4) -------------------------------------------------
indice = (1:N);
dx_50 = cos(phi_50).*d_50;
dy_50 = sin(phi_50).*d_50;

dx_100 = cos(phi_100).*d_100;
dy_100 = sin(phi_100).*d_100;

figure
subplot(2,2,1)
scatter(indice, dx_50, 1)
title('Dx 50')
subplot(2,2,2)
scatter(indice, dy_50, 1)
title('Dy 50')
subplot(2,2,3)
scatter(indice, dx_100, 1)
title('Dx 100')
subplot(2,2,4)
scatter(indice, dy_100, 1)
title('Dy 100')

figure
subplot(2,1,1)
scatter(dx_50, dy_50, 1)
axis([0 60 1 30])
subplot(2,1,2)
scatter(dx_100, dy_100, 1)
axis([0 105 1 80])

%----------------------------------------------------

%5) -------------------------------------------------

%----------------------------------------------------

%6) -------------------------------------------------

%----------------------------------------------------

%7) -------------------------------------------------

%----------------------------------------------------

%8) -------------------------------------------------

%----------------------------------------------------

%9) -------------------------------------------------

%----------------------------------------------------
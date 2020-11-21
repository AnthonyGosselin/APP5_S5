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

phi0_15 = 15;
phi0_30 = 30;

% 1)-------------------------------------------------
% a)------------------------
N = 10000;
theta = rand(1, N) * 2 * pi;
% --> Faire le graphique des nombres générés?

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

% Calculs en degrés, puis convertit en radians:
% D0 = 50
d_50 = d0_50 + r_50.*cos(theta);
phi_50_1 = (phi0_50 + r_50.*sin(theta))/ 180 * pi;
phi_50_2 = (phi0_30 + r_50.*sin(theta))/ 180 * pi;

% D0 = 100
d_100 = d0_100 + r_100.*cos(theta);
phi_100_2 = (phi0_100 + r_100.*sin(theta))/ 180 * pi;
phi_100_1 = (phi0_15 + r_100.*sin(theta))/ 180 * pi;

%----------------------------------------------------

%4) -------------------------------------------------
indice = (1:N);
dx_50_1 = cos(phi_50_1).*d_50;
dy_50_1 = sin(phi_50_1).*d_50;

dx_50_2 = cos(phi_50_2).*d_50;
dy_50_2 = sin(phi_50_2).*d_50;


dx_100_1 = cos(phi_100_1).*d_100;
dy_100_1 = sin(phi_100_1).*d_100;

dx_100_2 = cos(phi_100_2).*d_100;
dy_100_2 = sin(phi_100_2).*d_100;

% Graphiques 1D ---
% 15 deg
figure
subplot(2,2,1)
scatter(indice, dx_50_1, 1)
title('Dx 50 (15 deg)')
subplot(2,2,2)
scatter(indice, dy_50_1, 1)
title('Dy 50 (15 deg)')
subplot(2,2,3)
scatter(indice, dx_100_1, 1)
title('Dx 100 (15 deg)')
subplot(2,2,4)
scatter(indice, dy_100_1, 1)
title('Dy 100 (15 deg)')
% 30 deg
figure
subplot(2,2,1)
scatter(indice, dx_50_2, 1)
title('Dx 50 (30 deg)')
subplot(2,2,2)
scatter(indice, dy_50_2, 1)
title('Dy 50 (30 deg)')
subplot(2,2,3)
scatter(indice, dx_100_2, 1)
title('Dx 100 (30 deg)')
subplot(2,2,4)
scatter(indice, dy_100_2, 1)
title('Dy 100 (30 deg)')

% Graphiques 2D ---
% 15 deg
figure
subplot(2,2,1)
scatter(dx_50_1, dy_50_1, 1)
title('Dx 50 (15 deg)')
axis([0 60 1 35])
subplot(2,2,2)
scatter(dx_100_1, dy_100_1, 1)
title('Dx 100 (15 deg)')
axis([0 120 1 80])
%30 deg
subplot(2,2,3)
scatter(dx_50_2, dy_50_2, 1)
title('Dy 50 (30 deg)')
axis([0 60 1 35])
subplot(2,2,4)
scatter(dx_100_2, dy_100_2, 1)
title('Dy 100 (30 deg)')
axis([0 120 1 80])

%----------------------------------------------------

%5) -------------------------------------------------
%%
funcs_to_plot = {dx_50_1, dx_50_2, dy_50_1, dy_50_2, dx_100_1, dx_100_2, dy_100_1, dy_100_2};
funcs_titles = {'Dx = 50 (15 deg)', 'Dx = 50 (30 deg)', 'Dy = 50 (15 deg)', 'Dy = 50 (30 deg)', 'Dx = 100 (15 deg)', 'Dx = 100 (30 deg)', 'Dy = 100 (15 deg)', 'Dy = 100 (30 deg)'};
             
%figure
funcs_moyenne = {}
funcs_ecart_type = {}
for i = [1:length(funcs_to_plot)]
    figure(99)
    subplot(2, 4, i);
    %subplot(1, 2, 1);
    data = cell2mat(funcs_to_plot(i));
    histogram(data);
    plot_title = funcs_titles(i);
    title(plot_title);
    axis([-5 120 0 900]);
    
    moyenne = mean(data);
    ecart_type = sqrt(var(data));
    
    funcs_moyenne(i) = mat2cell(moyenne, 1);
    funcs_ecart_type(i) = mat2cell(ecart_type, 1);
    
    % Freqence relative
    figure(98)
    subplot(2, 4, i);
    %subplot(1, 2, 2)
    [freq_abs, edges] = histcounts(data);
    bin_mdpt=(edges(2:end)+edges(1:(end-1)))/2;
    freq_rel = freq_abs/N;
    stem(bin_mdpt, freq_rel);
    title(strcat(plot_title, ' freq rel(%)'));
    axis([-5 120 0 0.09]);
    txt = {['Moyenne = ', num2str(moyenne)],['s = ', num2str(ecart_type)]};
    text(10,0.08,txt)
    
end
%----------------------------------------------------

%6) -------------------------------------------------
% AH BEN CALINE, WHAT A SURPRISE, ON DIRAIT DES COURBES NORMALES

% Les distributions Dx et Dy ressemblent tous à des Normales (avec moyennes
% et ecart type différents.

% Les valeurs des moyennes et ecarts-types estimés dans les arrays: 'funcs_moyenne'
% et 'funcs_ecart_type'

% (On pourrait plot les courbes théoriques par dessus les courbes
% d'échantillon: avec normpdf() )

%----------------------------------------------------

%7) -------------------------------------------------

%----------------------------------------------------

%8) -------------------------------------------------

%----------------------------------------------------

%9) -------------------------------------------------

%----------------------------------------------------
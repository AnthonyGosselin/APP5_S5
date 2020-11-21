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
% D0 = 50
d_50 = sqrt(d0_50^2 + 2*d0_50*r_50.*cos(theta) + r_50.^2);
phi_50_1 = phi0_15 + atand(r_50.*sin(theta)./(d0_50 + r_50.*cos(theta)));
phi_50_2 = phi0_30 + atand(r_50.*sin(theta)./(d0_50 + r_50.*cos(theta)));

% D0 = 100
d_100 = sqrt(d0_100^2 + 2*d0_100*r_100.*cos(theta) + r_100.^2);
phi_100_1 = phi0_15 + atand(r_100.*sin(theta)./(d0_100 + r_100.*cos(theta)));
phi_100_2 = phi0_30 + atand(r_100.*sin(theta)./(d0_100 + r_100.*cos(theta)));

% % Calculs en degrés, puis convertit en radians:
% % D0 = 50
% d_50 = d0_50 + r_50.*cos(theta);
% phi_50_1 = (phi0_50 + r_50.*sin(theta))/ 180 * pi;
% phi_50_2 = (phi0_30 + r_50.*sin(theta))/ 180 * pi;
% 
% % D0 = 100
% d_100 = d0_100 + r_100.*cos(theta);
% phi_100_2 = (phi0_100 + r_100.*sin(theta))/ 180 * pi;
% phi_100_1 = (phi0_15 + r_100.*sin(theta))/ 180 * pi;

%----------------------------------------------------

%4) -------------------------------------------------
% Valeurs Dx et Dy totales
indice = (1:N);
dx_50_1 = cosd(phi_50_1).*d_50;
dy_50_1 = sind(phi_50_1).*d_50;

dx_50_2 = cosd(phi_50_2).*d_50;
dy_50_2 = sind(phi_50_2).*d_50;

dx_100_1 = cosd(phi_100_1).*d_100;
dy_100_1 = sind(phi_100_1).*d_100;

dx_100_2 = cosd(phi_100_2).*d_100;
dy_100_2 = sind(phi_100_2).*d_100;


% Valeurs Dx et Dy déterministes
d0x_50_1 = cosd(phi0_15)*d0_50;
d0y_50_1 = sind(phi0_15)*d0_50;

d0x_50_2 = cosd(phi0_30)*d0_50;
d0y_50_2 = sind(phi0_30)*d0_50;

d0x_100_1 = cosd(phi0_15)*d0_100;
d0y_100_1 = sind(phi0_15)*d0_100;

d0x_100_2 = cosd(phi0_30)*d0_100;
d0y_100_2 = sind(phi0_30)*d0_100;


%Valeurs Dx et Dy aléatoires
dx_50_1_alea = dx_50_1 - d0x_50_1;
dy_50_1_alea = dy_50_1 - d0y_50_1;

dx_50_2_alea = dx_50_2 - d0x_50_2;
dy_50_2_alea = dy_50_2 - d0y_50_2;

dx_100_1_alea = dx_100_1 - d0x_100_1;
dy_100_1_alea = dy_100_1 - d0y_100_1;

dx_100_2_alea = dx_100_2 - d0x_100_2;
dy_100_2_alea = dy_100_2 - d0y_100_2;


% Plot des fonctions
%funcs_to_plot = {dx_50_1, dx_50_2, dy_50_1, dy_50_2, dx_100_1, dx_100_2, dy_100_1, dy_100_2};
funcs_to_plot = {dx_50_1_alea, dx_50_2_alea, dy_50_1_alea, dy_50_2_alea, ...
                 dx_100_1_alea, dx_100_2_alea, dy_100_1_alea, dy_100_2_alea};
funcs_titles = {'Dx = 50 (15 deg)', 'Dx = 50 (30 deg)', 'Dy = 50 (15 deg)', 'Dy = 50 (30 deg)', 'Dx = 100 (15 deg)', 'Dx = 100 (30 deg)', 'Dy = 100 (15 deg)', 'Dy = 100 (30 deg)'};

figure(123)
for i = [1:length(funcs_to_plot)]
    subplot(2, 4, i);
    data = cell2mat(funcs_to_plot(i));
    scatter(indice, data, 1);
    plot_title = funcs_titles(i);
    title(plot_title);
end

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

% 30 deg
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
% et ecart type différents).

% (On pourrait plot les courbes théoriques par dessus les courbes
% d'échantillon: avec normpdf() )

figure
for i = [1:length(funcs_to_plot)] 
    subplot(2, 4, i)
    data = cell2mat(funcs_to_plot(i));
    [freq_abs, edges] = histcounts(data);
    bin_mdpt=(edges(2:end)+edges(1:(end-1)))/2;
    freq_rel = freq_abs/N;
    stem(bin_mdpt, freq_rel);

    hold on
    mu = mean(data);
    sigma = sqrt(var(data));
    norm_X = (bin_mdpt(1):0.01:bin_mdpt(end));
    plot(norm_X, normpdf(norm_X, mu, sigma));
    
    plot_title = funcs_titles(i);
    title(plot_title);
end

% Hmm, ok la courbe normale ne semble pas correspondre à toutes les
% courbes... (mais juste par un facteur d'échelle tout le temps)
%----------------------------------------------------

%7) -------------------------------------------------

%----------------------------------------------------

%8) -------------------------------------------------
% Calcul de covariance, à la main
mat_50_1_prime = [dx_50_1_alea - mean(dx_50_1_alea); dy_50_1_alea - mean(dy_50_1_alea)];
mat_50_2_prime = [dx_50_2_alea - mean(dx_50_2_alea); dy_50_2_alea - mean(dy_50_2_alea)];
mat_100_1_prime = [dx_100_1_alea - mean(dx_100_1_alea); dy_100_1_alea - mean(dy_100_1_alea)];
mat_100_2_prime = [dx_100_2_alea - mean(dx_100_2_alea); dy_100_2_alea - mean(dy_100_2_alea)];

mat_cov_50_1 = mat_50_1_prime * mat_50_1_prime' / N;
mat_cov_50_2 = mat_50_2_prime * mat_50_2_prime' / N;
mat_cov_100_1 = mat_100_1_prime * mat_100_1_prime' / N;
mat_cov_100_2 = mat_100_2_prime * mat_100_2_prime' / N;

% Covariance avec fonction matlab
mat_cov_50_1_fct = cov(dx_50_1_alea, dy_50_1_alea);
mat_cov_50_2_fct = cov(dx_50_2_alea, dy_50_2_alea);
mat_cov_100_1_fct = cov(dx_100_1_alea, dy_100_1_alea);
mat_cov_100_2_fct = cov(dx_100_2_alea, dy_100_2_alea);

%----------------------------------------------------

%9) -------------------------------------------------
%


%----------------------------------------------------
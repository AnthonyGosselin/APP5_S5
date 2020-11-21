%%LABORATOIRE
%% Numero 1
close all
clear all
clc

for N = [100 1000 10000]

    random_5to5 = rand(1, N)*10 - 5;

    edges = (-5:1:5);
    figure
    histogram(random_5to5, edges);

    a = -5;
    b = 5;

    mu = (a+b)/2;
    sigma = sqrt(((b-a)^2)/12);

    xbar = mean(random_5to5);
    s = sqrt(var(random_5to5));
    
    mean_err = (mu - xbar)^2;
    var_err = (sigma-s)^2;
    
    figure
end

%% Numero 2
clear cll
close all
clc

for N = [100 1000 10000]

    distribution = randn(1, N)*2 + 10;

    figure
    histogram(distribution);

    mu = 10;
    sigma = 2;

    xbar = mean(distribution);
    s = sqrt(var(distribution));
    
    err_mean = (mu - xbar)^2;
    err_var = (sigma-s)^2;
    
end

%% Numero 3
clear all
close all
clc


for N = [100 1000 10000]

    U1 = rand(1,N);
    U2 = rand(1,N);
    
    mu = 10;
    sigma = 2;
    
    X1 = mu + sigma*cos(2*pi*U1).*sqrt(-2*log(U2));
    X2 = mu + sigma*sin(2*pi*U1).*sqrt(-2*log(U2));
    
    figure
    histogram(X1);
    figure
    histogram(X2);
    
    xbar1 = mean(X1);
    s1 = sqrt(var(X1));
    
    err_mean1 = (mu - xbar1)^2;
    err_var1 = (sigma-s1)^2;
    
    xbar2 = mean(X2);
    s2 = sqrt(var(X2));
    
    err_mean2 = (mu - xbar2)^2;
    err_var2 = (sigma-s2)^2;
end

%% Numero 4
clear all
close all
clc

for N = [100 1000 10000]
    p = (0:(1/N):1);
    x = sqrt(-log(1-(2*p-1).^2)/ sqrt(pi/8));

    index = 1;
    for xval = x
        if p(index) < 0.5
           x(index) = -xval; 
        end
        index = index+1;
    end

    figure
    plot(x, p)

    %b
    distribution = rand(1,N);

    cloche = [];
    index = 1;
    for val = distribution
       val = floor(val * N)+1;
       cloche(index) = x(val);
       index = index+1;
    end
    
    sigma = 2;
    mu = 10;
    new_cloche = cloche*sigma + mu;

    figure
    histogram(new_cloche)
    
    xbar = mean(new_cloche);
    s = sqrt(var(new_cloche));
    
    err_mean = (mu - xbar)^2;
    err_var = (sigma-s)^2;
end

%% Numero 5
clear all
close all
clc

for N = [100 1000 10000 100000]
    p = (0:(1/N):1);
    beta = 10;
    
    x = -beta*log(1-p);

    figure
    plot(x, p)

    %b
    distribution = rand(1,N);

    cloche = [];
    index = 1;
    for val = distribution
       val = floor(val * N)+1;
       cloche(index) = x(val);
       index = index+1;
    end
    figure
    histogram(cloche)
    
    xbar = mean(cloche);
    s = sqrt(var(cloche));
    
     err_mean = (beta - xbar)^2;
     err_var = (beta-s)^2;
end

%% Partie 2 - Cercle
clear all
close all
clc

N = 10000000;
x = rand(1,N);
y = rand(1,N);

d = sqrt(x.^2 + y.^2);

N2 = 0;
for dval = d
    if dval < 1
        N2 = N2 + 1;
    end
end

approx_pi = 4 * N2/N

%% Partie 2 - Sphere
close all
clear all
clc

N = 1000000;
x = rand(1,N);
y = rand(1,N);
z = rand(1,N);

d = sqrt(x.^2 + y.^2 + z.^2);

N3 = 0;
for dval = d
    if dval < 1
        N3 = N3 + 1;
    end
end

approx_pi = 8 * N3/N



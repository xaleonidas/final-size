% Initialisation:

n = 50; % Initial susceptibles
m = 1; % Initial infectives
lambda = 2; % Rate of infectious contact

%Laplace transform:
syms x;
I = 4*x*exp(-2*x); % Infectious time probability distribution
phi = laplace(I);

%res = double(subs(phi, lambda/n));

probs_2 = zeros(1,n+1);

% l=0:
probs_2(1) = double(subs(phi, lambda))^m;

% Other cases:
for l = 1:n
    total = 0;
    for k = 0:l-1
        res = nchoosek(n-k, l-k) * probs_2(k+1) * double(subs(phi, lambda*(n-l)/n))^(l-k);
        total = total + res;
    end
    probs_2(l+1) = nchoosek(n,l) * double(subs(phi, lambda*(n-l)/n))^(l+m) - total;
end

plot(probs_2);
title('Numeric Example - n=50,m=1,\lambda=2')
xlabel('Final size Z') 
ylabel('Probability') 
exportgraphics(gca,'selke_50_1_2.pdf','ContentType','vector')
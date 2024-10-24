% Initialisation:

n = 57; % Initial susceptibles
m = 2; % Initial infectives
lambda = 2; % Rate of infectious contact

% Make martix of probabilitities for each generation.
% s is susceptible number, i is infective number.

% rank 0:
prob_mat = zeros(n+1);
prob_mat(n+1, m+1) = 1;

for rank = 1:n+1
    old_prob_mat = prob_mat;
    prob_mat = zeros(n+1);
    for s = 0:n
        for i = 0:n
            if s+i <= n
                for j = 0:n
                    prob_mat(s+1,i+1) = prob_mat(s+1,i+1) + old_prob_mat(s+i+1, j+1) ...
                        * trans(i, s+i, j, lambda, n);
                end
            end
        end
    end
    rank
end

probs = flip(prob_mat(:,1));

function prob = trans(i, s, j, lambda, n)
    %Laplace transform:
    syms x;
    I = 4*x*exp(-2*x); % Infectious time probability distribution
    phi = laplace(I);
    temp_sum = 0;
    for k=0:i
        temp_sum = temp_sum + nchoosek(i, k) * (-1)^k * ...
            double(subs(phi, lambda*(s-i+k)/n))^j;
    end
    prob = temp_sum * nchoosek(s, i);
end
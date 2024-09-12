%% Max-Robust-Sum (MRS) Test Statistics Implementation
% The implementation includes the computation of the null distribution function through both simulation and analytical methods. 
% It also plots the results from both methods for empirical verification and theoretical analysis.

%% Numerical simulation
n = 20;  % Example value
m = 15;  % Example value
p = 2;   % Example value
j = m - p;

num_simulations = 20000;
stat = zeros(num_simulations, 1); 
for k = 1 : num_simulations
    data = exprnd(1, n, 1);
    sorted_data = sort(data, 'descend');
    stat_aux = sorted_data(j) / sum(sorted_data(m + 1 : end));
    stat(k) = stat_aux;
end
figure;
histogram(stat, 500, 'Normalization', 'pdf'); 
hold on;

%% Analytical derivation 
n = vpa(n);
m = vpa(m);

% Define the range for t > 1 / (n-m)
lb = double(1 / (n - m)) + 0.01;  % Example value
ub = 2;                           % Example value
t_values = linspace(lb, ub, 1000); 
t_values = vpa(t_values);

% Preallocate the array for results
MRS_pdf = vpa(zeros(size(t_values)));

% Calculate PDF
for i = 1 : length(t_values)
    % Standardization
    MRS_pdf(i) = MRS_df(t_values(i), n, m, p) / area(n, m, p);
end

% Plot the results
plot(t_values, MRS_pdf, 'k', 'LineWidth', 1.5);
xlabel('t');
ylabel('pdf');
title('Null distribution function of the MRS');
grid on;
xlim([0 ub]);

%%
function result = MRS_df(t, n, m, p)
    digits(50);  % Example value
    result = vpa(0);  % Initialization
    
    parfor r = 1 : p
        % Binomial coefficient for the first part
        binom_coeff1 = nchoosek(p - 1, r - 1);
        
        % First major term
        term1 = factorial(n - m) * ...
                ((1 + r * (1/(-m + n) - t) + m * t)^(-1 + m - n) - ...
                 (1 + m * t)^(-1 + m - n)) / r;

        % First sub-sum
        innerSum1 = vpa(0);
        for k = 1 : floor(1/t)
            binom_coeff2 = nchoosek(n - m, k);
            subSum1 = vpa(0);
            for l = 1 : n - m - 1
                subSum1 = subSum1 + (1 + l) * (r * (-(1/k) + t))^l * (1 + m * t)^(-2 - l);
            end
            innerSum1 = innerSum1 + (-1)^k * binom_coeff2 * k^(n - m - 1) * subSum1;
        end
        
        % Second sub-sum
        innerSum2 = vpa(0);
        for k = 1 : floor(1/t)
            binom_coeff2 = nchoosek(n - m, k);
            innerSum2 = innerSum2 + (-1)^k * binom_coeff2 * k^(n - m - 1) / (1 + m * t)^2;
        end
        
        % Third sub-sum
        innerSum3 = vpa(0);
        for k = (floor(1/t) + 1) : n - m - 1
            binom_coeff2 = nchoosek(n - m, k);
            innerSum3 = innerSum3 + (-1)^k * binom_coeff2 * k^(n - m - 1) / (1 + m * t + r * (1/k - t))^2;
        end
        
        % Fourth sub-sum
        innerSum4 = vpa(0);
        for k = 1 : n - m - 1
            binom_coeff2 = nchoosek(n - m, k);
            subSum2 = vpa(0);
            for l = 0 : n - m - 1
                subSum2 = subSum2 + (1 + l) * (r * (-(1/k) + 1/(-m + n)))^l * ...
                            (1 + r * (1/(-m + n) - t) + m * t)^(-2 - l);
            end
            innerSum4 = innerSum4 + (-1)^k * binom_coeff2 * k^(n - m - 1) * subSum2;
        end
        
        % Combine nested sums
        term2 = (-r)^(m - n) * factorial(n - m - 1) * (innerSum1 + innerSum2 + innerSum3 - innerSum4);
        
        % Add to the result
        result = result + (-1)^(p - r) * binom_coeff1 * (term1 + term2);
    end
end

%%
function result = area(n, m, p)
    digits(50);  % Example value
    result = vpa(0);  % Initialization
    
    parfor r = 1 : p
        % Binomial coefficient for the first part
        binom_coeff1 = nchoosek(p - 1, r - 1);
        
        % First major term
        term1 = factorial(n - m) * ...
                (((n/(-m + n))^(1 + m - n) / (n * (m - r))) - ...
                 ((-(n/(m - n)))^(1 + m - n) / (m * n))) / r;

        % First sub-sum
        innerSum1 = vpa(0);

        for k = 1 : n - m - 1
            binom_coeff2 = nchoosek(n - m, k);
            subSum1 = vpa(0);
            for l = 1 : n - m - 1
                subSum1 = subSum1 + (1 + l) * r^l * (-((k + m - n)^(1 + l) * (k * n)^(-l)) / ((1 + l) * (k + m) * n));
            end
            innerSum1 = innerSum1 + (-1)^k * binom_coeff2 * k^(n - m - 1) * subSum1;
        end
        
        % Second sub-sum
        innerSum2 = vpa(0);
        for k = 1 : n - m - 1
            binom_coeff2 = nchoosek(n - m, k);
            innerSum2 = innerSum2 + (-1)^k * binom_coeff2 * k^(n - m - 1) * (1/(k+m) - 1/n);
        end
        
        % Third sub-sum
        innerSum3 = vpa(0);
        for k = 1 : n - m - 1
            binom_coeff2 = nchoosek(n - m, k);
            innerSum3 = innerSum3 + (-1)^k * binom_coeff2 * k^(n - m - 1) * k / ((k+m) * (m-r));
        end
        
        % Fourth sub-sum
        innerSum4 = vpa(0);
        for k = 1 : n - m - 1
            binom_coeff2 = nchoosek(n - m, k);
            subSum2 = vpa(0);
            for l = 0:(n - m - 1)
                subSum2 = subSum2 + (1 + l) * (r * (-(1/k) + 1/(-m + n)))^l * ...
                            (n / (n-m))^(-1 - l) /  ((1+l) * (m-r));
            end
            innerSum4 = innerSum4 + (-1)^k * binom_coeff2 * k^(n - m - 1) * subSum2;
        end
        
        % Combine nested sums
        term2 = (-r)^(m - n) * factorial(n - m - 1) * (innerSum1 + innerSum2 + innerSum3 - innerSum4);
        
        % Add to the result
        result = result + (-1)^(p - r) * binom_coeff1 * (term1 + term2);
    end
end

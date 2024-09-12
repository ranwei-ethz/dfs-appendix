%% Sum-Robust-Sum (SRS) Test Statistics Implementation
% The implementation includes the computation of the null distribution function through both simulation and analytical methods. 
% It also plots the results from both methods for empirical verification and theoretical analysis.

%% Numerical simulation
n = 12;  % Example value
m = 6;   % Example value
r = 3;   % Example value
p = m - r - 1;  

num_simulations = 20000;
stats = zeros(num_simulations, 1); 
for k = 1:num_simulations
    data = exprnd(1, n, 1);
    sorted_data = sort(data, 'descend');
    stats_aux = sum(sorted_data(1 : r))/sum(sorted_data(m + 1 : end));
    stats(k) = stats_aux;
end
figure;
histogram(stats, 500, 'Normalization', 'pdf'); 
hold on;

%% Analytical derivation 
n = vpa(n);
m = vpa(m);

% Define the range for t > r / (n-m)
lb = double(r / (n - m)) + 0.01;  % Example value
ub = 10;                           % Example value
y_values = linspace(lb, ub, 1000); 
y_values = vpa(y_values);

% Preallocate the array for results
SRS_pdf = vpa(zeros(size(y_values)));

% Calculate the sum for each t
for i = 1:length(y_values)
    SRS_pdf(i) = SRS_df(y_values(i), n, m, r, p) / area(n, m, r, p);
end

% Plot the results
plot(y_values, SRS_pdf, 'k', 'LineWidth', 1.5);
xlabel('y');
ylabel('pdf');
title('Null distribution function of the SRS');
grid on;
xlim([0 10]);

%%
function result = SRS_df(y, n, m, r, p)
    digits(50);  % Set the precision to 50 decimal places
    result = vpa(0);  % Initialize the sum
    
    parfor i = 1:p
        
        % First inner sum (l) I
        innerSum11 = vpa(0);
        for l = 0:r-1
            innerSum11 = innerSum11 + (-((m - r) * (r / (m - n) + y)) / r)^l * ...
                        ((r + m * y - n * (1 + y)) / (m - n))^(-1 - l + m - n) * ...
                        factorial(n - m + l) / factorial(l);
        end

        % First term I
        term11 = (-(m - r))^(-r) * (innerSum11 - (1 + (m * y) / r)^(-1 + m - n) * factorial(n - m));

        % First inner sum (l) II
        innerSum12 = vpa(0);
        for l = 0:r-1
            innerSum12 = innerSum12 + (((i - m + r) * (r / (m - n) + y)) / r)^l * ...
                        ((r + m * y - n * (1 + y)) / (m - n))^(-1 - l + m - n) * ...
                        factorial(n - m + l) / factorial(l);
        end

        % First term II
        term12 = (i - m + r)^(-r) * ((1 + i / (n - m) + (m - i) * y / r)^(-1 + m - n) * factorial(n - m) - innerSum12);

        % First major term
        term1 = factorial(r - 1) / i * (term11 + term12); 

        % Manual split
        max_k = floor(r / y);

        % Second inner sum (j, h, l, k) I
        innerSum21 = vpa(0);
        for k = 1:max_k  

            % Sum over l
            innerSum21_l = vpa(0);
            for l = 1:n-m-1
                
                % Sum over h
                innerSum21_h = vpa(0);
                for h = 0:r-1
    
                    % Sum over j
                    innerSum21_j = vpa(0);
                    for j = 1:h+l
                        subTerm21_j1 = (-(k + m - n) * (m - r) / (k * (m - n)))^j * ...
                                    ((r + m * y - n * (1 + y)) / (m - n))^(-1 + h - j - r);
                        subTerm21_j2 = (((m - r) * (-r + k * y)) / (k * r))^j * ...
                                    (1 + (m * y) / r)^(-1 + h - j - r);
                        innerSum21_j = innerSum21_j + (factorial(r - h + j) / factorial(j)) * (subTerm21_j1 - subTerm21_j2);
                    end
    
                    % Combine terms
                    innerSum21_h = innerSum21_h + (-1)^h * nchoosek(r - 1, h) * (y/r - 1/k)^(r - 1 - h) * (m - r)^(-1 - h - l) * ...
                                   factorial(h + l) * (factorial(r - h) * (((r + m * y - n * (1 + y)) / (m - n))^(-1 + h - r) - ...
                                   (1 + (m * y)/r)^(-1 + h - r)) + innerSum21_j);
                end
                
                % Accumulate the result
                innerSum21_l = innerSum21_l + i^l / factorial(l) * innerSum21_h;
            end
    
            % Final accumulation into the outer result
            innerSum21 = innerSum21 + (-1)^k * nchoosek(n - m, k) * k^(n - m - 1) * innerSum21_l;
        end

        % Second inner sum (j, h, l, k) II
        innerSum22 = vpa(0);
        for k = max_k+1:n-m-1  

            % Sum over l
            innerSum22_l = vpa(0);
            for l = 1:n-m-1
                
                % Sum over h
                innerSum22_h = vpa(0);
                for h = 0:r-1
    
                    % Sum over j
                    innerSum22_j = vpa(0);
                    for j = 0:h+l
                        subTerm22_j = ((k + m - n) * (-m + r) / (k * (m - n)))^j * ...
                                   ((r + m * y - n * (1 + y)) / (m - n))^(-1 + h - j - r);
                        innerSum22_j = innerSum22_j + (factorial(r - h + j) / factorial(j)) * subTerm22_j;
                    end
    
                    % Combine terms
                    innerSum22_h = innerSum22_h + (-1)^h * nchoosek(r - 1, h) * (y/r - 1/k)^(r - 1 - h) * (m - r)^(-1 - h - l) * ...
                                   factorial(h + l) * (innerSum22_j - factorial(r - h) * ((k + m - r + k * y) / k)^(-1 + h - r));
                end
                
                % Accumulate the result
                innerSum22_l = innerSum22_l + i^l / factorial(l) * innerSum22_h;
            end
    
            % Final accumulation into the outer result
            innerSum22 = innerSum22 + (-1)^k * nchoosek(n - m, k) * k^(n - m - 1) * innerSum22_l;
        end

        % Second major term
        term2 = innerSum21 + innerSum22; 
     
        % Third inner sum (l, k) I
        innerSum31 = vpa(0);
        for k = 1:max_k

            innerSum31_l = vpa(0);
            for l = 0:r-1
                innerSum31_l = innerSum31_l + (1 + l) * (m - n)^2 * (((m - r) * (r + (m - n) * y)) / r)^l * (n - r - m * y + n * y)^(-2 - l);
            end
            innerSum31 = innerSum31 + (-1)^k * nchoosek(n - m, k) * k^(n - m - 1) * (-m + r)^(-r) * factorial(r - 1) * ...
                (r^2 / (r + m * y)^2 - innerSum31_l);
        end

        % Third inner sum (l, k) II
        innerSum32 = vpa(0);
        for k = max_k+1:n-m-1

            innerSum32_l = vpa(0);
            for l = 0:r-1
                innerSum32_l = innerSum32_l + ((1 + l) * k^2 * (((m - r) * (r - k * y)) / r)^l * (k + m - r + k * y)^(-2 - l) - ...
                    (1 + l) * (m - n)^2 * (((m - r) * (r + (m - n) * y)) / r)^l * (n - r - m * y + n * y)^(-2 - l));
            end
            innerSum32 = innerSum32 + (-1)^k * nchoosek(n - m, k) * k^(n - m - 1) * (-m + r)^(-r) * factorial(r - 1) * innerSum32_l;
        end

        term3 = innerSum31 + innerSum32;

        
        % Fourth inner sum (l, k)
        term4 = vpa(0);
        for k = max_k+1:n-m-1

            innerSum4_l = vpa(0);
            for l = 0:r-1
                innerSum4_l = innerSum4_l + (1 + l) * (-(((i - m + r) * (r - k * y)) / r))^l * (k + m - r + k * y)^(-2 - l);
            end
            term4 = term4 + (-1)^k * nchoosek(n - m, k) * k^(n - m + 1) * (i - m + r)^(-r) * factorial(r - 1) * ...
                (r^2 / (i * (r - k * y) + k * (r + m * y))^2 - innerSum4_l);
        end


        % Fifth inner sum (h, l, k)
        term5 = vpa(0);
        for k = 1:n-m-1

            innerSum5_l = vpa(0);
            for l = 0:n-m-1

                innerSum5_h = vpa(0);
                for h = 0:r-1
                    innerSum5_h = innerSum5_h + nchoosek(l + h + 1, h) * (-m + n)^(2 + l) * ...
                        (-((i - m + r) * (r + (m - n) * y) / r))^h * (n - r - m * y + n * y)^(-2 - h - l);
                end

                innerSum5_l = innerSum5_l + (l + 1) * (-i * (1/k + 1/(m - n)))^l * (i - m + r)^(-r) * factorial(r - 1) * ...
                    ((1 + i/(-m + n) + ((-i + m) * y)/r)^(-2 - l) - innerSum5_h);
            end

            term5 = term5 + (-1)^k * nchoosek(n - m, k) * k^(n - m - 1) * innerSum5_l;
        end

        
        % Final sum combination
        result = result + (-1)^(p - i) * nchoosek(p - 1, i - 1) * r^(-1 + r) * ...
            (term1 + (-i)^(m - n) * factorial(n - m - 1) * (term2 + term3 + term4 - term5));
    end
end

%%
function result = area(n, m, r, p)
    digits(50);  % Example value
    result = vpa(0);  % Initialization
    
    parfor i = 1:p
        
        % First inner sum (l) I
        innerSum11 = vpa(0);
        for l = 0:r-1
            innerSum11 = innerSum11 + (1 - m / n)^(n - m) * (1 - m / r)^l;
        end

        % First term I
        term11 = (-(m - r))^(-r) * (innerSum11 - (n / (n - m))^(1 + m - n) * r * (n - m) / (m * n));

        % First inner sum (l) II
        innerSum12 = vpa(0);
        for l = 0:r-1
            innerSum12 = innerSum12 + (1 - m / n)^(n - m) * ((i - m + r) / r)^l;
        end

        % First term II
        term12 = (i - m + r)^(-r) * ((n / (n - m))^(1 + m - n) * r * (n - m) / ((m - i) * n) - innerSum12);

        % First major term
        term1 = factorial(r - 1) / i * (term11 + term12); 

        % Second inner sum (j, h, l, k) I
        innerSum21 = vpa(0);
        for k = 1:n-m-1  

            % Sum over l
            innerSum21_l = vpa(0);
            for l = 1:n-m-1
                
                % Sum over h
                innerSum21_h = vpa(0);
                for h = 0:r-1
    
                    % Sum over j
                    innerSum21_j = vpa(0);
                    for j = 1:h+l

                        innerSum21_q = vpa(0);
                        for q = 0:r-h-1
                            innerSum21_q = innerSum21_q + nchoosek(r - h - 1, q) * (-r / k - (n - r) / (n - m))^q * ...
                                ((r / k + (n - r) / (n - m))^(-1 - j - q) - (n / (n - m))^(-1 - j - q)) / (-1 - j - q);
                        end

                        subTerm21_j1 = (-(k + m - n) * (m - r) / (k * (m - n)))^j * r^(1 + h - r);
                        subTerm21_j2 = ((k * n) / (k + m - n))^(h - j - r) * (m -r)^j * k * r / (k + m) / (h - j - r);
                        innerSum21_j = innerSum21_j + (factorial(r - h + j) / factorial(j)) * (subTerm21_j1 * innerSum21_q - subTerm21_j2);
                    end
    
                    % Combine terms
                    innerSum21_h = innerSum21_h + (-1)^h * nchoosek(r - 1, h) * (m - r)^(-1 - h - l) * factorial(h + l) * ...
                                   (factorial(r - h - 1) * k * r * ((k + m - n) / (k * n))^(r - h) * ...
                                   (1 / (k + m) - (n - m) / (k * n - (k + m - n) * r)) + innerSum21_j);
                end
                
                % Accumulate the result
                innerSum21_l = innerSum21_l + i^l / factorial(l) * innerSum21_h;
            end
    
            % Final accumulation into the outer result
            innerSum21 = innerSum21 + (-1)^k * nchoosek(n - m, k) * k^(n - m - 1) * innerSum21_l;
        end

        % Second inner sum (j, h, l, k) II
        innerSum22 = vpa(0);
        for k = 1:n-m-1  

            % Sum over l
            innerSum22_l = vpa(0);
            for l = 1:n-m-1
                
                % Sum over h
                innerSum22_h = vpa(0);
                for h = 0:r-1
    
                    % Sum over j
                    innerSum22_j = vpa(0);
                    for j = 0:h+l
                        innerSum22_j = innerSum22_j + (((k + m - n) * (m - r)) / (k * (n - r) + r * (n - m)))^j;
                    end
    
                    % Combine terms
                    innerSum22_h = innerSum22_h + (-1)^h * nchoosek(r - 1, h) * (m - r)^(-1 - h - l) * factorial(h + l) * k * ...
                                   factorial(r - h - 1) * r^(1 + h - r) * (-1 / (k + m) + (n - m) / (k * n - (k + m - n) * r) * innerSum22_j);
                end
                
                % Accumulate the result
                innerSum22_l = innerSum22_l + i^l / factorial(l) * innerSum22_h;
            end
    
            % Final accumulation into the outer result
            innerSum22 = innerSum22 + (-1)^k * nchoosek(n - m, k) * k^(n - m - 1) * innerSum22_l;
        end

        % Second major term
        term2 = innerSum21 + innerSum22; 
     
        % Third inner sum (l, k) I
        innerSum31 = vpa(0);
        for k = 1:n-m-1  

            innerSum31_l = vpa(0);
            for l = 0:r-1
                innerSum31_l = innerSum31_l + (n - m) * (m - r)^l * ((k + m - n) / (k * n - (k + m - n) * r))^(l + 1);
            end
            innerSum31 = innerSum31 + (-1)^k * nchoosek(n - m, k) * k^(n - m - 1) * (-m + r)^(-r) * factorial(r) / n * ...
                (innerSum31_l - (k + m - n) / (k + m));
        end

        % Third inner sum (l, k) II
        innerSum32 = vpa(0);
        for k = 1:n-m-1

            innerSum32_l = vpa(0);
            for l = 0:r-1
                innerSum32_l = innerSum32_l + (1 - m / r)^l * ((n - m) / n * ((1 - (k * n) / (k * n - (k + m - n) * r))^(1 + l) - 1) + k / (k + m));
            end
            innerSum32 = innerSum32 + (-1)^k * nchoosek(n - m, k) * k^(n - m - 1) * (-m + r)^(-r) * factorial(r - 1) * innerSum32_l;
        end

        term3 = innerSum31 + innerSum32;

        
        % Fourth inner sum (l, k)
        term4 = vpa(0);
        for k = 1:n-m-1

            innerSum4_l = vpa(0);
            for l = 0:r-1
                innerSum4_l = innerSum4_l + ((i - m + r) / r)^l;
            end
            term4 = term4 + (-1)^k * nchoosek(n - m, k) * k^(n - m) * (i - m + r)^(-r) * factorial(r - 1) / (k + m) * ...
                (-r / (i - m) - innerSum4_l);
        end


        % Fifth inner sum (h, l, k)
        term5 = vpa(0);
        for k = 1:n-m-1

            innerSum5_l = vpa(0);
            for l = 0:n-m-1

                innerSum5_h = vpa(0);
                for h = 0:r-1
                    innerSum5_h = innerSum5_h + nchoosek(l + h + 1, h) * ((-m + n) / n)^(1 + l) * ...
                        ((i - m + r) / r)^h * factorial(h) * factorial(l) / factorial(h + l + 1);
                end

                innerSum5_l = innerSum5_l + (l + 1) * (-i * (1/k + 1/(m - n)))^l * (i - m + r)^(-r) * factorial(r - 1) * ...
                    (((1 - m / n)^l * (m - n) * r) / ((1 + l) * (i - m) * n) - innerSum5_h);
            end

            term5 = term5 + (-1)^k * nchoosek(n - m, k) * k^(n - m - 1) * innerSum5_l;
        end

        
        % Final sum combination
        result = result + (-1)^(p - i) * nchoosek(p - 1, i - 1) * r^(-1 + r) * factorial(n - m - 1) *...
            (term1 + (-i)^(m - n) * (term2 + term3 + term4 - term5));
    end
end

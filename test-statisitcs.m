%% Simulation and Density Estimation of the Test Statistic
% This script simulates and computes the density of a self-normalized CUSUM-based test statistic.
% Additionally, it calculates the p-value for observed realizations of the test statistics
% in an upper-tailed hypothesis test.

%% Number of Monte Carlo iterations
num_it = 4000;
Data = zeros(1, num_it); % Preallocate for efficiency

for l = 1:num_it
    
    %% Simulating Brownian Motion
    N = 2000; % Number of steps
    dt = 1 / N; % Time step size
    dW = sqrt(dt) * randn(1, N); % Wiener process increments
    BrownianMotion = cumsum([0, dW]); % Construct Brownian motion path
    
    %% Compute Self-Normalized CUSUM Statistic
    B = zeros(1, N + 1);
    D = zeros(1, N + 1);
    Num = zeros(1, N + 1);
    
    for k = 1:N+1
        % Compute numerator: deviation from scaled Brownian motion
        Num(k) = BrownianMotion(k) - (k * dt) * BrownianMotion(N+1);
        
        % Compute the integral-based denominator
        for s = 1:k
            B(s) = BrownianMotion(s) - (s/k) * BrownianMotion(k);
        end
        for s = k+1:N+1
            B(s) = BrownianMotion(s) - BrownianMotion(k) - ((s-k)/(N-k)) * (BrownianMotion(N+1) - BrownianMotion(k));
        end
        B = B.^2;
        
        % Compute integrals for denominator
        integral1 = dt * (B(1) + B(k)) * 0.5 + dt * sum(B(2:k-1));
        integral2 = dt * (B(k) + B(N+1)) * 0.5 + dt * sum(B(k+2:N));
        D(k) = sqrt(integral1 + integral2);
    end
    
    % Compute test statistic for current iteration
    a = abs(Num) ./ D;
    Data(l) = max(a);
end

%% Step 1: Plot histogram of the simulated test statistic
histogram(Data, 30)
title('Histogram of Simulated Test Statistic')
xlabel('Test Statistic Value')
ylabel('Frequency')

%% Step 2: Estimate Probability Density Function (PDF)
[g, u] = ksdensity(Data); % Kernel density estimation

%% Step 3: Plot the estimated PDF
figure;
plot(u, g, 'LineWidth', 2);
xlabel('X');
ylabel('Density');
hold on;

% Plot mean as a vertical reference line
xline(mean(Data), 'r--', 'LineWidth', 1); % Red dashed line at mean value

title('Estimated Probability Density Function');
legend('$\sup_{\tau \in [0,1]} \frac{|B(\tau)-\tau B(1)|}{ \left[\int_0^\tau \left( B(s) - \frac{s}{ \tau} B(\tau) \right)^2 \ ds+ \int_{\tau}^1 \left( B(s) - B(\tau) -\frac{s-\tau}{1-\tau} \left( B(1)-B(\tau) \right)\right)^2 \ ds\right]^{1/2}}\;.$', 'Interpreter', 'latex', 'FontSize', 17);

%% Step 4: Compute quantiles for significance levels
p = [0.975, 0.95, 0.025]; % Commonly used quantiles
Quantiles = quantile(Data, p);

%% Step 5: Compute p-value for an upper-tailed test
observed_value = 1.5469; % Observed test statistic (replace with actual value if needed)
p_value = mean(Data >= observed_value); % P-value computation

%% Display results
fprintf('Observed value: %.4f\n', observed_value);
fprintf('Upper-tailed p-value: %.6f\n', p_value);
% This script simulates the theoretical distribution of a self-normalized statistic
% using Brownian motion paths. 
% The statistic is calculated as the maximum of a self-normalized Brownian bridge
% divided by a normalization factor involving integrals of squared deviations.
%
% Key outputs:
% - 'Data': Simulated values of the statistic for num_it Monte Carlo iterations.
% - 'Data1' and 'Data2': Values of the test statistics for two Moving Average models.
%
% Notes:
% - The statistic is defined as:
%   max |B(t) - tB(1)| / sqrt([integrals over deviations of the Brownian bridge])
% - %innovations can be used to modify the innovation process if needed.
% - The script also plots estimated PDFs and histograms for the simulated data.

%% Simulation Parameters
num_it = 4000; % Number of Monte Carlo simulations
Data = zeros(1, num_it); % Array to store results for each iteration

%% Simulate the Statistic for num_it Iterations
for l = 1:num_it
    N = 2000;           % Number of steps in the Brownian motion
    dt = 1 / N;         % Time step
    dW = sqrt(dt) * randn(1, N); % Generate N increments for the Brownian motion
    BrownianMotion = cumsum([0, dW]); % Construct the Brownian motion path

    % Initialize arrays for computations
    B = zeros(1, N + 1);
    D = zeros(1, N + 1);
    Num = zeros(1, N + 1);

    % Compute the self-normalized statistic
    for k = 1:N+1
        % Calculate the numerator
        Num(k) = BrownianMotion(k) - (k * dt) * BrownianMotion(N+1);

        % Calculate the Brownian bridge deviations
        for s = 1:k
            B(s) = BrownianMotion(s) - (s/k) * BrownianMotion(k);
        end
        for s = k+1:N+1
            B(s) = BrownianMotion(s) - BrownianMotion(k) - ...
                   ((s-k)/(N-k)) * (BrownianMotion(N+1) - BrownianMotion(k));
        end
        B = B .^ 2; % Square the deviations

        % Calculate the integrals
        integral1 = dt * (B(1) + B(k)) * 0.5 + dt * sum(B(2:k-1));
        integral2 = dt * (B(k) + B(N+1)) * 0.5 + dt * sum(B(k+2:N));
        D(k) = sqrt(integral1 + integral2); % Normalization factor
    end

    % Calculate the statistic for this iteration
    a = abs(Num) ./ D;
    Data(l) = max(a);
end

%% Values of the Test Statistics for Moving Average Models
% test_MA_correct calculates the test statistic for testing
% X_t = eps_t + a * eps_{t-1}, with Moving Average models.

a = 0.4;
b = 0.4;
w = 0.5; % Proportion of time series used for testing

[Data1, ~] = test_MA_correct(a, b, 1000, w);
[Data2, ~] = test_MA_correct(a, b + 0.3, 1000, w);

%% Plotting the Results

% Kernel density estimation for the simulated data
[g, u] = ksdensity(Data);

% Figure 1: PDF and histogram for Data1
figure(1);
plot(u, g, 'LineWidth', 2);
hold on;
histogram(Data1, 'Normalization', 'pdf');
xlabel('SC_{n_b}', 'FontSize', 16);
ylabel('Density', 'FontSize', 16);
title('PDF and Histogram for Data1', 'FontSize', 16);
grid on;

% Figure 2: PDF and histogram for Data2
figure(2);
plot(u, g, 'LineWidth', 2);
hold on;
histogram(Data2, 'Normalization', 'pdf');
xlabel('SC_{n_b}', 'FontSize', 16);
ylabel('Density', 'FontSize', 16);
title('PDF and Histogram for Data2', 'FontSize', 16);
grid on;

% This script creates N stationary AR(1) time series, each containing T data points.
% A structural change occurs at time 'tau', where the AR(1) coefficient changes:
% Xt = phi * X_{t-1} + eps_t,   for t = 1, ..., tau
% Xt = phi2 * X_{t-1} + eps_t,  for t = tau+1, ..., T
% 
% Key variables:
% - T: Number of points in each time series.
% - eps_t: Innovations (errors) whose distribution can change (e.g., Gaussian, Laplace, etc.).
% - tau: Location of the Change Point, expressed as a fraction in [0,1].
% - phi, phi2: AR(1) coefficients before and after the change point.
% - %innovations: This command allows the user to modify the distribution of the innovations.
%
% Outputs:
% - 'data': Self-normalized CUSUM values for each time series.
% - 'v': Frequency (in percentage) of rejecting the null hypothesis ("no change")
%        for different change point locations and heights.
%
% The results provide insight into the power of detecting changes in AR(1) processes.

%% Initialization
delta = 0.001;         % Block size control parameter (must be in (0, 0.5))
T = 500;               % Number of points in each time series
m = floor(T^(0.5 + delta)); % Block size
N = 1000;              % Number of simulations
x = zeros(N, T);       % Initialize AR(1) process
v = zeros(6, 3);       % Store results for different heights and change locations

count1 = 0;            % Counter for different heights (h values)

%% Loop over different heights (h) and change point locations (tau)
for h = [0, 0.1, 0.2, 0.3, 0.4, 0.5]
    count1 = count1 + 1;
    phi = 0.4;            % AR(1) coefficient before the change
    phi2 = phi + h;       % AR(1) coefficient after the change
    count2 = 0;           % Counter for different tau values
    
    for tau = [0.1, 0.25, 0.5]  % Change point locations
        count2 = count2 + 1;
        T1 = floor(tau * T);    % Calculate index of the change point
        
        % Loop over each simulation
        for i = 1:N
            % Generate AR(1) process before and after the change point
            for t = 2:T1
                % %innovations: Customize the distribution of innovations here
                % Example options:
                % innovation = trnd(2);         % t-distributed innovations (df = 2)
                % innovation = randn;          % Gaussian (0,1)
                u = rand() - 0.5;              % Uniform random number in [-0.5, 0.5]
                innovation = -4 * sign(u) * log(1 - 2 * abs(u)); % Laplace(0, 4)
                x(i, t) = phi * x(i, t-1) + innovation;
            end
            for t = (T1+1):T
                % %innovations: Customize the distribution of innovations here
                % innovation = trnd(2);         % t-distributed innovations (df = 2)
                % innovation = randn;          % Gaussian (0,1)
                u = rand() - 0.5;
                innovation = -4 * sign(u) * log(1 - 2 * abs(u)); % Laplace(0, 4)
                x(i, t) = phi2 * x(i, t-1) + innovation;
            end
        end

        %% CUSUM Calculation
        data = zeros(1, N); % Store CUSUM values for all simulations
        for i = 1:N
            [alpha, ~] = turning_rate_byepoch(m, cumsum(x(i, :)));
            n = length(alpha);
            M = mean(alpha);
            Num = zeros(1, n);
            R = zeros(1, n);
            D = zeros(1, n);
            for k = 1:n
                Num(k) = abs(sum(alpha(1:k)) - k * M);
                for r = 1:k
                    R(r) = sum(alpha(1:r)) - (r/k) * sum(alpha(1:k));
                end
                for r = k+1:n
                    R(r) = sum(alpha(k+1:r)) - ((r-k)/(n-k)) * sum(alpha(k+1:n));
                end
                D(k) = sqrt(mean(R.^2));
            end
            data(i) = max(Num ./ D); % Store the maximum self-normalized value
        end

        %% Test Null Hypothesis (H0: "No Change")
        rej = zeros(1, N); % Number of rejections
        for i = 1:N
            if data(i) >= 6.3325
                rej(i) = 1;
            end
        end
        v(count1, count2) = mean(rej) * 100; % Percentage of rejections
    end
end

%% Plot Results
figure;
plot([0, 0.1, 0.2, 0.3, 0.4, 0.5], v(:, 1), '-o', 'LineWidth', 2, 'Color', [0 0.4470 0.7410], ...
    "MarkerFaceColor", [0.9290 0.6940 0.1250]); % First tau
hold on;
plot([0, 0.1, 0.2, 0.3, 0.4, 0.5], v(:, 2), '-o', 'LineWidth', 2, 'Color', [0.4660 0.6740 0.1880], ...
    "MarkerFaceColor", [0.9290 0.6940 0.1250]); % Second tau
plot([0, 0.1, 0.2, 0.3, 0.4, 0.5], v(:, 3), '-o', 'LineWidth', 2, 'Color', [0.6350 0.0780 0.1840], ...
    "MarkerFaceColor", [0.9290 0.6940 0.1250]); % Third tau

xlabel('h = \phi_2 - \phi_1', 'FontSize', 18);
ylabel('Frequency of Detected Changes (%)', 'FontSize', 16);
legend('\tau = 0.1', '\tau = 0.25', '\tau = 0.5', 'Location', 'Best', 'FontSize', 20);
title(sprintf('Data points n = %d', T), 'FontSize', 16);
grid on;

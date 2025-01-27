function [data, v] = test_MA_correct(a, b, N, w)
% test_MA_correct generates N stationary Gaussian time series using an MA(1) scheme 
% with a structural change at time tau and calculates the self-normalized CUSUM statistic.
%
% Time series structure:
% Xt = eps_t + a * eps_{t-1}, for t = 1, ..., tau
% Xt = eps_t + b * eps_{t-1}, for t = tau+1, ..., T
%
% Key inputs:
% - a, b: MA(1) coefficients before and after the change point, respectively.
% - N: Number of time series realizations.
% - w: Fractional location of the change point (tau = floor(w * T)).
%
% Outputs:
% - data: Self-normalized CUSUM values for each time series.
% - v: Percentage of rejections of the null hypothesis ("no change").
%
% Method:
% - Simulates N time series with a structural change at tau.
% - Computes self-normalized CUSUM statistics for turning rates.
% - Evaluates the null hypothesis ("no change").
%
% Notes:
% - eps_t are Gaussian innovations (stationary process).
% - The rejection threshold for the null hypothesis is hardcoded as 6.4425.

%% Initialization
delta = 0.001; % Controls block size (must be in (0, 0.5))
T = 5000;      % Number of points in each time series
m = floor(T^(0.5 + delta)); % Block size
tau = floor(T * w); % Change point location (in indices)

% Create the covariance structure for MA(1) process with a change point
e1 = ones(tau, 1);     % Pre-change coefficients
e2 = ones(T - tau, 1); % Post-change coefficients
S = [spdiags([a * e1, e1], -1:0, tau, tau), sparse(tau, T - tau); ...
     sparse(T - tau, tau), spdiags([b * e2, e2], -1:0, T - tau, T - tau)];
S = [sparse(T, 1), S];
S(1, 1) = a;
S(tau + 1, tau) = b;

% Sample N time series from a multivariate normal distribution
meanVector = zeros(T, 1);
sample = mvnrnd(meanVector, S * S', N);

%% Compute Self-Normalized CUSUM Statistics
data = zeros(1, N); % Initialize array to store statistics
for i = 1:N
    [alpha, ~] = turning_rate_byepoch(m, sample(i, :));
    n = length(alpha);
    M = mean(alpha);
    Num = zeros(1, n);
    R = zeros(1, n);
    D = zeros(1, n);

    % Calculate CUSUM numerator and normalization factor
    for k = 1:n
        Num(k) = abs(sum(alpha(1:k)) - k * M); % Numerator
        for r = 1:k
            R(r) = sum(alpha(1:r)) - (r/k) * sum(alpha(1:k));
        end
        for r = k+1:n
            R(r) = sum(alpha(k+1:r)) - ((r-k)/(n-k)) * sum(alpha(k+1:n));
        end
        D(k) = sqrt(mean(R.^2)); % Denominator
    end

    % Store the maximum CUSUM statistic
    data(i) = max(Num ./ D);
end

%% Test Null Hypothesis (H0: "no change")
rej = zeros(1, N); % Array to count rejections
threshold = 6.4425; % Predefined rejection threshold
for i = 1:N
    if data(i) >= threshold
        rej(i) = 1;
    end
end

% Calculate percentage of rejections
v = mean(rej) * 100;
fprintf('Number of rejections: %.2f%%\n', v);

end

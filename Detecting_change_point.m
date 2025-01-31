function [alpha, time] = Detecting_change_point(X)
% Detecting Change Point in a Time Series
% This function applies a CUSUM-based method to detect change points in a time series.
% The x-axis of the plots is customized according to the specific example in the paper.
% In the paper, the 5th patient from the Terzano dataset is used, and the segment 13,200,000 to 14,400,000 is analyzed.

%% Step 1: Compute CUSUM-based Change Point Detection
% X is the time series for which we want to detect the change point.

delta = 0.001; % Parameter controlling block size, must be in (0, 1/2)
m = floor((length(X))^(0.5 + delta)); % Compute block size
[alpha, ~] = turning_rate_byepoch(m, cumsum(X)); % Compute turning rate

n = length(alpha);
M = mean(alpha);
Num = zeros(1, n);
R = zeros(1, n);
D = zeros(1, n);

for k = 1:n
    Num(k) = abs(sum(alpha(1:k)) - k * M);
    
    for r = 1:k
        R(r) = sum(alpha(1:r)) - (r / k) * sum(alpha(1:k));
    end
    for r = k+1:n
        R(r) = sum(alpha(k+1:r)) - ((r-k) / (n-k)) * sum(alpha(k+1:n));
    end
    
    D(k) = sqrt(mean(R.^2));
end

[data, time] = max(Num ./ D);

disp('The index is: ')
time = time * m;
disp(time)
plot(Num ./ D)

%% Step 2: Test Hypothesis H0 (No Change)
if data >= 6.3325
    disp('Reject H0: Change detected')
end
disp('Test statistic value: ')
disp(data)

%% Step 3: Plot Results
figure
subplot(2, 1, 1)
plot(X, 'Color', 'k')
xlabel('TIME', 'FontSize', 16);
ylabel('\muV', 'FontSize', 18);
ylim([-30, 30])
xlim([1, length(X)])

% Convert detected time to actual time reference
result = time_for_observation(13200000 + time, '17:58:48');
set(gca, 'xtick', [1 time length(X)], 'xticklabel', {'01:08:2'; result; '01:47:33'})
hold on

% Mark detected change point
xline(time, 'r--', 'LineWidth', 3);
title('Transition from REM to S2', 'FontSize', 20);

subplot(2, 1, 2)
plot(alpha, 'Color', 'k')
hold on
plot(mean(alpha(1:time/m)) * ones(time/m), 'b--', 'LineWidth', 2)
hold on
plot(time/m:n, mean(alpha(time/m:end)) * ones(n+1-time/m), 'm--', 'LineWidth', 2)
xlabel('n_b', 'FontSize', 16);
ylabel('q', 'FontSize', 18);
ylim([0, 0.07])
xlim([1, length(alpha)])

set(gca, 'xtick', [1 time/m length(alpha)], 'xticklabel', {1; time/m; n})
hold on

% Mark detected change point in turning rate plot
xline(time, 'r--', 'LineWidth', 3);
title('Turning Rate', 'FontSize', 20);
end

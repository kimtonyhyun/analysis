function [thresh, stats] = estimate_baseline_threshold(trace)

stats = compute_trace_stats(trace);

tr_lower = trace(trace <= stats.mode);
sigma = std(tr_lower - stats.mode);

% Convert standard deviation of the half-normal distribution 
% to that of the normal distribution
sigma = sigma / (1 - 2/pi);
thresh = stats.mode + 3*sigma;

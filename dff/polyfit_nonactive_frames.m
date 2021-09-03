function [baseline, info] = polyfit_nonactive_frames(trace, threshold, padding, order)
% Use 'parse_active_frames' to identify periods of activity in the calcium
% trace. Perform 'polyfit' to the periods of non-activity, in order to
% estimate the trace baseline.
%

if isempty(threshold)
    m = min(trace);
    M = max(trace);
    threshold = m + 1/2*(M-m);
end

num_frames = length(trace);
t = linspace(0, 1, num_frames);

active_frames = parse_active_frames(trace > threshold, padding);
active_frames = frame_segments_to_list(active_frames);

nonactive_frames = true(num_frames, 1);
nonactive_frames(active_frames) = false;

if (sum(nonactive_frames) < order)
    cprintf('red', 'Insufficient number of nonactive frames for fitting!\n');
    
    % As fallback, fit to all frames
    nonactive_frames = true(num_frames, 1);
    threshold = max(trace);
end

x = t(nonactive_frames);
y = trace(nonactive_frames);

p = polyfit(x, y, order);
baseline = polyval(p, t);

info.threshold = threshold;
info.padding = padding;
info.order = order;
info.nonactive_frames = nonactive_frames;
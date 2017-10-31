function draw_constant_path_analysis(ds, cell_idx)
% Used by 'view_raster_touch'

trials = ds.get_switch_trials;

pre_trials = trials.constant_pre;
post_trials = trials.constant_post;

x_pre = 1:length(pre_trials);
x_post = x_pre(end) + (1:length(post_trials));
x_range = [1 x_post(end)];

pre = compute_features(ds, cell_idx, pre_trials);
post = compute_features(ds, cell_idx, post_trials);


% Behavioral trial times
subplot(4,2,2);
stem(x_pre, pre.times, 'b');
hold on;
stem(x_post, post.times, 'r');
hold off;
xlim(x_range);
ylabel('Trial times (s)');
ylim(compute_ylim(pre.times, post.times));
grid on;
title('Constant path: Pre (blue) vs. Post (red)');

% Fluorescence averages
subplot(4,2,4);
stem(x_pre, pre.mean_fluorescence, 'b');
hold on;
stem(x_post, post.mean_fluorescence, 'r');
hold off;
xlim(x_range);
ylabel('Avg fluorescence');
ylim(compute_ylim(pre.mean_fluorescence, post.mean_fluorescence));
grid on;

% Event amplitude sum
subplot(4,2,6);
stem(x_pre, pre.event_sum, 'b');
hold on;
stem(x_post, post.event_sum, 'r');
hold off;
xlim(x_range);
ylabel('\Sigma Event amplitudes');
ylim(compute_ylim(pre.event_sum, post.event_sum));
grid on;

subplot(4,2,8);
stem(x_pre, pre.num_events, 'b');
hold on;
stem(x_post, post.num_events, 'r');
hold off;
xlabel('Trial index');
xlim(x_range);
ylabel('Event counts');
ylim(compute_ylim(pre.num_events, post.num_events));
grid on;

end % draw_constant_path_analysis

function y_range = compute_ylim(pre_vals, post_vals)
    all_vals = [pre_vals(:); post_vals(:)];
    m = min(min(all_vals), 0);
    M = max(all_vals);
    y_range = [m M] + 0.1*(M-m)*[-1 1];
end
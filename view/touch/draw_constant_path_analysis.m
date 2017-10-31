function draw_constant_path_analysis(ds, cell_idx)
% Used by 'view_raster_touch'

trials = ds.get_switch_trials;

pre_trials = trials.constant_pre;
post_trials = trials.constant_post;

pre = compute_features(ds, cell_idx, pre_trials);
post = compute_features(ds, cell_idx, post_trials);


% Behavioral trial times
subplot(4,2,2);
stem(pre_trials, pre.times, 'b');
hold on;
stem(post_trials, post.times, 'r');
hold off;
% xlabel('Trial index');
xlim([1 ds.num_trials]);
ylabel('Trial times (s)');
grid on;
title('Constant path: Pre (blue) vs. Post (red)');

% Fluorescence averages
subplot(4,2,4);
stem(pre_trials, pre.mean_fluorescence, 'b');
hold on;
stem(post_trials, post.mean_fluorescence, 'r');
hold off;
% xlabel('Trial index');
xlim([1 ds.num_trials]);
grid on;
ylabel('Avg fluorescence');

% Event amplitude sum
subplot(4,2,6);
stem(pre_trials, pre.event_sum, 'b');
hold on;
stem(post_trials, post.event_sum, 'r');
hold off;
% xlabel('Trial index');
xlim([1 ds.num_trials]);
ylabel('\Sigma Event amplitudes');
grid on;

subplot(4,2,8);
stem(pre_trials, pre.num_events, 'b');
hold on;
stem(post_trials, post.num_events, 'r');
hold off;
xlabel('Trial index');
xlim([1 ds.num_trials]);
ylabel('Event counts');
grid on;
function draw_constant_path_analysis(ds, cell_idx)
% Used by 'view_raster_touch'

trials = ds.get_switch_trials;

pre_trials = trials.constant_pre;
post_trials = trials.constant_post;
num_pre = length(pre_trials);
num_post = length(post_trials);

% Behavioral trial times
pre_times = [ds.trials(pre_trials).time];
post_times = [ds.trials(post_trials).time];
subplot(4,2,2);
stem(pre_trials, pre_times, 'b');
hold on;
stem(post_trials, post_times, 'r');
hold off;
% xlabel('Trial index');
xlim([1 ds.num_trials]);
ylabel('Trial times (s)');
grid on;
title('Constant path: Pre (blue) vs. Post (red)');

% Comparison of fluorescence averages
pre_avg_fluorescence = zeros(1,num_pre);
post_avg_fluorescence = zeros(1,num_post);
for k = 1:num_pre
    trial_idx = pre_trials(k);
    tr = ds.get_trace(cell_idx, trial_idx);
    pre_avg_fluorescence(k) = mean(tr);
end
for k = 1:num_post
    trial_idx = post_trials(k);
    tr = ds.get_trace(cell_idx, trial_idx);
    post_avg_fluorescence(k) = mean(tr);
end

subplot(4,2,4);
stem(pre_trials, pre_avg_fluorescence, 'b');
hold on;
stem(post_trials, post_avg_fluorescence, 'r');
hold off;
all_fluorescence = [pre_avg_fluorescence, post_avg_fluorescence];
min_f = min(all_fluorescence);
max_f = max(all_fluorescence);
f_range = [min_f max_f] + 0.1*(max_f - min_f)*[-1 1];
% xlabel('Trial index');
xlim([1 ds.num_trials]);
ylim(f_range);
grid on;
ylabel('Avg fluorescence');

% Sum of event amplitudes
pre_event_sum = zeros(1,num_pre);
post_event_sum = zeros(1,num_post);

pre_event_count = zeros(1,num_pre);
post_event_count = zeros(1,num_post);

for k = 1:num_pre
    trial_idx = pre_trials(k);
    eventdata = ds.get_events(cell_idx, trial_idx);
    if ~isempty(eventdata)
        pre_event_sum(k) = sum(eventdata(:,3));
        pre_event_count(k) = size(eventdata,1);
    end
end
for k = 1:num_post
    trial_idx = post_trials(k);
    eventdata = ds.get_events(cell_idx, trial_idx);
    if ~isempty(eventdata)
        post_event_sum(k) = sum(eventdata(:,3));
        post_event_count(k) = size(eventdata,1);
    end
end

subplot(4,2,6);
stem(pre_trials, pre_event_sum, 'b');
hold on;
stem(post_trials, post_event_sum, 'r');
hold off;
% xlabel('Trial index');
xlim([1 ds.num_trials]);
ylabel('\Sigma Event amplitudes');
grid on;

subplot(4,2,8);
stem(pre_trials, pre_event_count, 'b');
hold on;
stem(post_trials, post_event_count, 'r');
hold off;
xlabel('Trial index');
xlim([1 ds.num_trials]);
ylabel('Event counts');
grid on;
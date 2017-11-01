function draw_constant_path_analysis(ds, cell_idx)
% Used by 'view_raster_touch'

trials = ds.get_switch_trials;

pre_trials = trials.constant_pre;
post_trials = trials.constant_post;

true_trial_inds = [pre_trials; post_trials];

x_pre = 1:length(pre_trials);
x_post = x_pre(end) + (1:length(post_trials));
x_range = [1 x_post(end)];

pre = compute_features(ds, cell_idx, pre_trials);
post = compute_features(ds, cell_idx, post_trials);

% Behavioral trial times
subplot(4,2,2);
draw_stem(pre.times, post.times);
ylabel('Trial times (s)');
title('Constant path: Pre (blue) vs. Post (red)');

% Fluorescence averages
subplot(4,2,4);
draw_stem(pre.mean_fluorescence, post.mean_fluorescence);
ylabel('Mean fluorescence');

% Event amplitude sum
subplot(4,2,6);
draw_stem(pre.event_sum, post.event_sum);
ylabel('\Sigma Event amplitudes');

% Event counts
subplot(4,2,8);
draw_stem(pre.num_events, post.num_events);
xlabel('Trial index');
ylabel('Event counts');

    function draw_stem(pre_vals, post_vals)
        stem(x_pre, pre_vals, 'b.', 'ButtonDownFcn', @select_trial);
        hold on;
        stem(x_post, post_vals, 'r.', 'ButtonDownFcn', @select_trial);
        hold off;
        xlim(x_range);
        tick_inds = x_range(1):5:x_range(end);
        xticks(tick_inds);
        xticklabels(num2cell(true_trial_inds(tick_inds)));
        ylim(compute_ylim(pre_vals, post_vals));
        grid on;
    end % draw_stem

    function select_trial(~, e)
        x = e.IntersectionPoint(1);
        t = true_trial_inds(x);
        msgbox(sprintf('Trial %d', t));
    end

end % draw_constant_path_analysis

function y_range = compute_ylim(pre_vals, post_vals)
    all_vals = [pre_vals(:); post_vals(:)];
    m = min(min(all_vals), 0);
    M = max(all_vals);
    
    if (m ~= M)
        y_range = [m M] + 0.1*(M-m)*[-1 1];
    else
        y_range = m + [-0.5 0.5];
    end
end
function draw_path_subrasters(ds, cell_idx, raster_scale)
% Used by 'view_raster_touch'

sd = ds.switchdata;
trials = ds.get_switch_trials;

constant_path_title = sprintf('Constant path (%s)', sd.constant_path_start);
changing_path_title = sprintf('Changing path (%s)', sd.changing_path_start);

% First column is for CONSTANT path
%------------------------------------------------------------
subplot(5,4,[3 7]);
% Note that 'aligned_time' is common to all traces
[constant_pre_traces, aligned_time] = draw_subraster(trials.constant_pre, 'PRE-switch');
title(constant_path_title);

subplot(5,4,[11 15]);
constant_post_traces = draw_subraster(trials.constant_post, 'POST-switch');
title(constant_path_title);

subplot(5,4,19);
draw_pre_post_averages(constant_pre_traces, constant_post_traces);
title(constant_path_title);

% Second column is for CHANGING path
%------------------------------------------------------------
subplot(5,4,[4 8]);
changing_pre_traces = draw_subraster(trials.changing_pre, 'PRE-switch');
title(changing_path_title);

subplot(5,4,[12 16]);
changing_post_traces = draw_subraster(trials.changing_post, 'POST-switch');
title(changing_path_title);

subplot(5,4,20);
draw_pre_post_averages(changing_pre_traces, changing_post_traces);
title(changing_path_title);

    function [raster, aligned_time] = draw_subraster(trial_inds, y_label)
        % Wrapper around 'plot_cell_raster' to apply some common formatting
        % for the current context
        [raster, aligned_time] = ds.plot_cell_raster(cell_idx, 'inds', trial_inds);
        set(gca, 'CLim', raster_scale);
        N = length(trial_inds);
        yticks(1:N);
        yticklabels(num2cell(trial_inds));
        ylabel(sprintf('%s (%d)', y_label, N));
        xlabel('');
    end % draw_subraster

    function draw_pre_post_averages(pre_traces, post_traces)
        % Compute the mean and the SEM
        pre_mean = mean(pre_traces,1);
        pre_err = std(pre_traces,[],1)/sqrt(size(pre_traces,1));
        post_mean = mean(post_traces,1);
        post_err = std(post_traces,[],1)/sqrt(size(post_traces,1));
        
        shadedErrorBar(aligned_time, pre_mean, pre_err, 'b');
        hold on;
        shadedErrorBar(aligned_time, post_mean, post_err, 'r');
        hold off;
        xlim(aligned_time([1 end]));
        ylim(raster_scale);
        grid on;
        xlabel('Frames relative to gate close');
        ylabel('Fluorescence');
    end % draw_pre_post_averages

end % draw_path_subrasters
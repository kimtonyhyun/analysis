function plot_trial_times(ds)

trial_indices = double(ds.trial_indices);
num_frames_per_trial = trial_indices(:,4) - trial_indices(:,1) + 1;
full_trial_times = num_frames_per_trial / ds.trace_fps;

% "Running time" is the period between onset of motion and gate close
movement_onset_frames = [ds.trials.movement_onset_frame]';
running_frames = (trial_indices(:,3) - trial_indices(:,1)) -...
                  movement_onset_frames + 1;
running_times = running_frames / ds.trace_fps;

ax1 = subplot(211);
plot(full_trial_times, '.-');
grid on;
xrange = [0.5 ds.num_trials+0.5];
xlim(xrange);
ylim(compute_yrange(full_trial_times));
ylabel('Full trial time [s]');

ax2 = subplot(212);
plot(running_times, 'r.-');
grid on;
xlim(xrange);
xlabel('Trial');
ylabel('Running time [s]');
yrange = compute_yrange(running_times);

% Draw trial correctness over the running time plot
corr_height = 0.5;
for k = 1:ds.num_trials
    if ds.trials(k).correct
        corr_color = 'g';
    else
        corr_color = 'r';
    end
    rectangle('Position', [k-0.5 yrange(2) 1 corr_height],...
              'FaceColor', corr_color);
end
ylim([yrange(1) yrange(2)+corr_height]);

linkaxes([ax1, ax2], 'x');

subplot(211); % Give focus to the topmost panel

end % plot_trial_times

function y_range = compute_yrange(times)
    y_min = min(times);
    y_max = max(times);
    y_range = [y_min y_max] + 0.1*(y_max-y_min)*[-1 1];
end
function view_trials_in_trace_by_color(time, trace, frame_indices)
% Displays the ICA trace color-coded by trial
% 2015 01 31 Tony Hyun Kim

num_trials = size(frame_indices,1);

colors = 'kbr';
for trial_idx = 1:num_trials
    trial_frames = frame_indices(trial_idx,1):...
                   frame_indices(trial_idx,end);
    plot(time(trial_frames),trace(trial_frames),...
         colors(mod(trial_idx,length(colors))+1));
    hold on;
end

xlim([time(1) time(end)]);
trace_delta = max(trace) - min(trace);
ylim([min(trace) max(trace)] + 0.1*trace_delta*[-1 1]);
xlabel('Time [s]');
ylabel('Signal [a.u.]');
grid on;
hold off;
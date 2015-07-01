function view_superimposed_trials_in_trace(trace, frame_indices)
% Displays the ICA trace as a superposition of trials

num_trials = size(frame_indices,1);

colors = 'kbr';
for trial_idx = 1:num_trials
    trial_frames = frame_indices(trial_idx,1):...
                   frame_indices(trial_idx,end);
    tr = trace(trial_frames);
    
    % Map the time onto the interval [0 1]
    plot(linspace(0,1,length(tr)), tr,...
         colors(mod(trial_idx,length(colors))+1));
    hold on;
end
xlabel('Phase of trial');
trace_delta = max(trace) - min(trace);
ylim([min(trace) max(trace)] + 0.1*trace_delta*[-1 1]);
ylabel('Signal [a.u.]');
grid on;
hold off;
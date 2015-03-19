function view_trace(time, trace, frame_indices)
% Displays the ICA trace color-coded by trial and as superposition

clf;
subplot(3,1,[2 3]);
view_superimposed_trials_in_trace(trace, frame_indices);

% Display the trace
subplot(3,1,1); % Ending on 311 makes it so that subsequent calls to plot
                %   e.g. 'title' go to the top panel
view_trials_in_trace_by_color(time, trace, frame_indices);

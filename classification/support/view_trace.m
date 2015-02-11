function view_trace(time, trace, frame_indices)
% Displays the ICA trace color-coded by trial and as superposition
% 2015 01 31 Tony Hyun Kim

% Params
mad_threshold_scale = 8;

% Trace properties
mad = compute_mad(trace);
thresh = mad_threshold_scale*mad;

clf;
subplot(3,1,[2 3]);
% set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]); % Maximize figure
view_superimposed_trials_in_trace(time, trace, frame_indices);
hold on;
plot([0 1], thresh*[1 1], 'r--', 'LineWidth', 2);
hold off;

% Display the trace
subplot(3,1,1); % Ending on 311 makes it so that subsequent calls to plot
                %   e.g. 'title' go to the top panel
view_trials_in_trace_by_color(time, trace, frame_indices);
hold on;
plot([time(1) time(end)], thresh*[1 1], 'r--', 'LineWidth', 2);
hold off;
function view_detailed_trial_touch(ds, cell_idx, trial_idx, varargin)

h_fig = [];
for i = 1:length(varargin)
    vararg = varargin{i};
    if ischar(vararg)
        switch lower(vararg)
            case 'fig'
                h_fig = varargin{i+1};
        end
    end
end

if isempty(h_fig)
    h_fig = figure;
else
    figure(h_fig);
    clf;
end

% Trial-specific data
%------------------------------------------------------------
[trial_frames, alignment_frame] = get_frames_relative_to_close(ds, trial_idx);
Mb = ds.get_behavior_trial(trial_idx); % Behavior movie
trace = ds.trials(trial_idx).traces(cell_idx,:);
num_frames_in_trial = length(trace);

% Image of cell
%------------------------------------------------------------
subplot(3,2,1);
imagesc(ds.cells(cell_idx).im);
axis image;
title(sprintf('Cell %d (%s)', cell_idx, ds.cells(cell_idx).label));
colormap jet; freezeColors;

% Zoom into the selected trial in the full raster
%------------------------------------------------------------
trial_window = 10;
trial_window_start = max(1, trial_idx-trial_window);
trial_window_end = min(ds.num_trials, trial_idx+trial_window);
gui.raster = subplot(3,2,[3 5]);
ds.plot_cell_raster(cell_idx, 'draw_correct');
scale = get(gca, 'CLim');
ylim([trial_window_start-0.5 trial_window_end+0.5]);
x_ends = get(gca, 'XLim');
line(x_ends, (trial_idx-0.5)*[1 1], 'Color', 'w', 'LineWidth', 2);
line(x_ends, (trial_idx+0.5)*[1 1], 'Color', 'w', 'LineWidth', 2);
title('All trials');
colormap jet; freezeColors;

set(gui.raster, 'ButtonDownFcn', @jump_to_trial);

% Show trace
%------------------------------------------------------------
gui.trace = subplot(3,2,2);
plot(trial_frames(1):trial_frames(4), trace, 'LineWidth', 2, 'HitTest', 'off');
xlim(trial_frames([1 4]));
ylim(scale);
grid on;
title(sprintf('Trial %d', trial_idx));
xlabel('Frames relative to gate close');
ylabel('Signal [a.u.]');
hold on;
plot(trial_frames(2)*[1 1], scale, 'k--', 'HitTest', 'off'); % Open-gate
plot(trial_frames(3)*[1 1], scale, 'k--', 'HitTest', 'off'); % Close-gate
gui.trace_bar = plot(trial_frames(1)*[1 1], scale, 'r', 'HitTest', 'off'); % Vertical bar
gui.trace_dot = plot(trial_frames(1), trace(1), 'r.',... % Dot
         'MarkerSize', 24,...
         'HitTest', 'off');

set(gui.trace, 'ButtonDownFcn', @update_frame);
     
% Show behavior movie
%------------------------------------------------------------
gui.behavior = subplot(3,2,[4 6]);
gui.behavior_image = imagesc(Mb(:,:,1));
set(gca, 'XTick', []);
set(gca, 'YTick', []);
axis image;
colormap gray;
title(sprintf('Frame 1 of %d', num_frames_in_trial));

% If tracking data loaded, overlay the positional information
if (ds.is_tracking_loaded)
    hold on;
    centroids = ds.trials(trial_idx).centroids;
    plot(centroids(:,1), centroids(:,2), '.-');
    gui.behavior_pos = plot(centroids(1,1), centroids(1,2), 'ro');
end

    function update_frame(~, e)
        sel_frame = round(e.IntersectionPoint(1));
        
        sel_frame = sel_frame + alignment_frame;
        sel_frame = max(sel_frame, 1);
        sel_frame = min(sel_frame, num_frames_in_trial);
        
        % Update visuals
        set(gui.trace_bar, 'XData', (sel_frame-alignment_frame)*[1 1]);
        set(gui.trace_dot, 'XData', (sel_frame-alignment_frame), 'YData', trace(sel_frame));
        set(gui.behavior_image, 'CData', Mb(:,:,sel_frame));
        
        subplot(gui.behavior);
        title(sprintf('Frame %d of %d', sel_frame, num_frames_in_trial));
        
        if (ds.is_tracking_loaded)
            centroid = ds.trials(trial_idx).centroids(sel_frame, :);
            set(gui.behavior_pos, 'XData', centroid(1), 'YData', centroid(2));
        end
    end % update_frame

    function jump_to_trial(~, e)
        t = round(e.IntersectionPoint(2));
        t = max(1, t);
        t = min(t, ds.num_trials);
        
        if (t ~= trial_idx)
            view_detailed_trial_touch(ds, cell_idx, t, 'fig', h_fig);
        end
    end % jump_to_trial

% Navigation controls
%------------------------------------------------------------
back_btn = uicontrol('Style', 'pushbutton',...
    'String', '<<',...
    'Units', 'normalized',...
    'Position', [0 0 0.05 1],...
    'Callback', @back_to_raster);

prev_trial_btn = uicontrol('Style', 'pushbutton',...
    'String', '^^',...
    'Units', 'normalized',...
    'Position', [0.95 0.5 0.05 0.5],...
    'Callback', {@trial_button_clicked, trial_idx-1});
if trial_idx == 1
    prev_trial_btn.Enable = 'off';
end

next_trial_btn = uicontrol('Style', 'pushbutton',...
    'String', 'vv',...
    'Units', 'normalized',...
    'Position', [0.95 0.0 0.05 0.5],...
    'Callback', {@trial_button_clicked, trial_idx+1});
if trial_idx == ds.num_trials
    next_trial_btn.Enable = 'off';
end

    function back_to_raster(~, ~)
        view_raster_touch(ds, cell_idx, 'fig', h_fig);
    end

    function trial_button_clicked(~, ~, trial_idx)
        view_detailed_trial_touch(ds, cell_idx, trial_idx, 'fig', h_fig);
    end

end % view_detailed_trial_touch

function [trial_frames, alignment_frame] = get_frames_relative_to_close(ds, trial_idx)
    trial_frames = double(ds.trial_indices(trial_idx, :)); % [start open-gate close-gate end]
    trial_frames = trial_frames - trial_frames(1) + 1;
    alignment_frame = trial_frames(3); % Align to close of gate

    trial_frames = trial_frames - alignment_frame;
end % get_frames_relative_to_close


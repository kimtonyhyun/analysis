function events = detect_events(ds, cell_idx, fps, varargin)
% Detect calcium events in the fluorescence trace interactively. Core
% computations are:
%   - Determine default thresholds (via 'estimate_baseline_sigma.m'),
%   - Compute events (via 'find_events_in_trials.m').
%
% The application allows for real-time visualization of detected events as
% a function of the algorithmic parameters.

cutoff_freq = 4; % Hz. From Wagner et al. 2019

use_prompt = true;
hfig = [];
M = [];
movie_clim = [];

for j = 1:length(varargin)
    vararg = varargin{j};
    if ischar(vararg)
        switch lower(vararg)
            case 'h' % Figure handle
                hfig = varargin{j+1};
            case 'noprompt'
                use_prompt = false;
            case 'cutoff'
                cutoff_freq = varargin{j+1};
            case {'m', 'movie'}
                M = varargin{j+1};
            case {'clim', 'movie_clim'}
                % To avoid recomputing movie_clim each time, which is a
                % costly operation.
                movie_clim = varargin{j+1};
        end
    end
end

% See if DaySummary already has events computed for this cell, or compute
% events from scratch.
%------------------------------------------------------------
trace_orig = ds.get_trace(cell_idx);

events = ds.cells(cell_idx).events;
if ~isempty(events)
    fprintf('  Using existing event detection results from DaySummary\n');
    [trace, stats] = preprocess_trace(fps, events.info.cutoff_freq);
else
    % No event data found in DaySummary
    events = struct('info', [], 'data', []);
    events.info = struct('num_frames', length(trace_orig),...
                         'cutoff_freq', [],...
                         'baseline', [],...
                         'sigma', [],...
                         'threshold', [],...
                         'merge_threshold', 0.05,...
                         'amp_threshold', 0.1);

    [trace, stats] = preprocess_trace(fps, cutoff_freq);
    events.info.threshold = events.info.baseline + 5*events.info.sigma;
    compute_events();
end

init_info = events.info;

% Set up GUI
%------------------------------------------------------------
if ~use_prompt
    return;
end

if ~isempty(M) && isempty(movie_clim)
    movie_clim = compute_movie_scale(M);
end

if isempty(hfig)
    hfig = figure;
end
gui = setup_gui(hfig, stats, trace_orig);

% GUI state
state.x_anchor = 1;
state.x_range = min(20*fps, gui.num_frames); % By default, show 20 s windows
state.show_orig = true;
state.show_dots = false;
state.show_trials = (ds.num_trials > 1);
state.show_neighbors = false;
state.sel_event = 0;
state.last_requested_trial = 0;
update_gui_state(gui, state);

if state.show_trials
    show_trial(1, gui);
end

% Interaction loop
%------------------------------------------------------------
prompt = 'Detector >> ';
while (use_prompt)
    resp = strtrim(input(prompt, 's'));
    val = str2double(resp);

    if (~isnan(val) && isreal(val)) % Is a number
        if (1 <= val) && (val <= ds.num_trials)
            show_trial(val, gui);
        end
    else % Not a number
        resp = lower(resp);
        if isempty(resp)
            resp = 'n';
        end
        switch (resp(1))
            case 'q' % "quit"
                ds_events = ds.cells(cell_idx).events;
                if isempty(ds_events)
                    ds.cells(cell_idx).events = events;
                else
                    if ~isequal(ds_events.info, events.info)
                        resp2 = strtrim(input('  Overwrite existing event results in DaySummary? (y/N) >> ', 's'));
                        switch lower(resp2)
                            case 'y'
                                ds.cells(cell_idx).events = events;
                        end
                    end
                end
                break;
                
            case {'s', 'w'} % Save event detection parameters and results to DaySummary
                ds.cells(cell_idx).events = events;
                fprintf('  Events saved to DaySummary\n');

            % Visualization
            %------------------------------------------------------------
            case {'i', 'z'} % zoom in
                x_center = state.x_anchor + 1/2 * state.x_range;
                state.x_range = 0.5*state.x_range;
                state.x_anchor = x_center - 1/2 * state.x_range;
                draw_local_window(gui, state);
                
            case 'o' % zoom out
                x_center = state.x_anchor + 1/2 * state.x_range;
                state.x_range = 2*state.x_range;
                state.x_anchor = x_center - 1/2 * state.x_range;
                draw_local_window(gui, state);
                
            case 'r' % toggle "raw"/original trace
                state.show_orig = ~state.show_orig;
                update_gui_state(gui, state);
                
            case 'd' % show dots
                state.show_dots = ~state.show_dots;
                update_gui_state(gui, state);
                
            case 'b' % show trials
                state.show_trials = ~state.show_trials;
                update_gui_state(gui, state);
                
            case 'x' % show neighbor correlations
                state.show_neighbors = ~state.show_neighbors;
                update_gui_state(gui, state);
                
            case 'n' % next trial
                if (ds.num_trials > 1) && (state.last_requested_trial < ds.num_trials)
                    show_trial(state.last_requested_trial+1, gui);
                end
                
            case 'm' % show even with minimum amplitude
                if ~isempty(events.data)
                    [~, ind] = min(events.data(:,3));
                    select_event(ind, gui);
                    
                    event_time = events.data(ind,2);
                    state.x_anchor = event_time - 1/2*state.x_range;
                    draw_local_window(gui, state);
                end
                
            % Modify processing parameters    
            %------------------------------------------------------------
            case 'e' % "Exclude": Increase the 'amp_threshold' to exclude currently selected event
                if state.sel_event == 0
                    fprintf('  Please first select an event to exclude\n');
                else
                    normamps = events.data(:,3) / max(events.data(:,3)); % amp_threshold is normalized
                    selected_amp = normamps(state.sel_event);
                    
                    normamps = sort(normamps, 'ascend');
                    ind = find(normamps>selected_amp, 1, 'first');
                    next_largest_amp = normamps(ind);
                    if isempty(next_largest_amp)
                        fprintf('  There are no events with a larger amplitude\n');
                    else
                        new_amp_threshold = 1/2*(selected_amp + next_largest_amp);
                        set_thresholds([], [], new_amp_threshold, gui);
                    end
                end
                
            case 'f' % Set lowpass filter cutoff frequency
                cf = str2double(resp(2:end));
                if cf > 0
                    fprintf('  Set cutoff frequency to %.1f\n', cf);
                    [trace, stats] = preprocess_trace(fps, cf);
                    
                    update_traces(gui);
                    draw_histogram(gui, stats);
                    set_thresholds([], [], [], gui); % Recomputes events
                else
                    fprintf('  Invalid cutoff frequency request "%s"\n', resp);
                end
                
            case 't' % reset parameters
                fprintf('  Reset event processing parameters to initial values\n');
                events.info = init_info;
                [trace, stats] = preprocess_trace(fps, events.info.cutoff_freq);
                update_traces(gui);
                draw_histogram(gui, stats);
                set_thresholds([], [], [], gui);

            otherwise
                fprintf('  Sorry, could not parse "%s"\n', resp);
        end
    end
end % Main interaction loop

    % Supplementary functions
    %------------------------------------------------------------
    function gui = setup_gui(hfig, stats, trace_orig)
        figure(hfig);
        clf(hfig);
        
        num_frames = length(trace_orig);
        trace_display_range = compute_display_range(trace_orig);
                
        % Display parameters kept around for convenience
        gui.hfig = hfig;
        gui.num_frames = num_frames;
        gui.trace_display_range = trace_display_range;
        
        % Setup the GLOBAL trace plot
        %------------------------------------------------------------
        gui.global = subplot(2,8,1:5);
        gui.global_rect = rectangle('Position',[-1 trace_display_range(1) 1 diff(trace_display_range)],...
                  'EdgeColor', 'none',...
                  'FaceColor', 'c', 'HitTest', 'off');
        hold on;
        gui.global_trace = plot(trace, 'k', 'HitTest', 'off');
        gui.global_thresh = plot([1 num_frames], -Inf*[1 1], 'm--', 'HitTest', 'off');
        gui.global_auto = plot(-Inf, -1, 'm.', 'HitTest', 'off');
        gui.global_manual = plot(-Inf, -1, 'r.', 'HitTest', 'off');
        hold off;
        box on;
        xlim([1 num_frames]);
        ylim(trace_display_range);
        xlabel('Frame');
        ylabel('Fluorescence (a.u.)');
        
        % Setup the HISTOGRAM plot
        %------------------------------------------------------------
        gui.histogram = subplot(2,8,6);
        gui.histogram_data = semilogy(1, 1, 'k.', 'HitTest', 'off');
        xlim(trace_display_range);
        hold on;
        gui.histogram_mode = plot(-Inf*[1 1], [0 1], '--', 'Color', 0.5*[1 1 1], 'HitTest', 'off');
        gui.histogram_thresh = plot(-Inf*[1 1], [0 1], 'm--', 'HitTest', 'off');
        hold off;
        ylabel('Trace histogram');
        view(-90, 90); % Rotate plot to match the global trace plot

        % Setup the post-trough height CDF
        %------------------------------------------------------------
        gui.post_amps = subplot(4,4,4);
        gui.post_cdf = plot(-1, -1, 'm.-', 'HitTest', 'off');
        hold on;
        gui.post_cdf_sel_event = plot(-1, -1, 'mo', 'HitTest', 'off');
        gui.post_cdf_amp_threshold = plot(-1*[1 1], [0 1], 'k', 'HitTest', 'off');
        hold off;
        xlim([0 0.5]);
        ylim([0 1]);
        xticks(0:0.1:0.5);
        yticks(0:0.1:1);
        grid on;
        xlabel('Post-amp / (max(trace) - baseline)');
        ylabel('Post-amp CDF');
        
        % Setup the pre-trough height CDF (i.e. event amplitude)
        %------------------------------------------------------------
        gui.pre_amps = subplot(4,4,8);
        gui.pre_cdf = plot(-1, -1, 'm.-', 'HitTest', 'off');
        hold on;
        gui.pre_cdf_sel_event = plot(-1, -1, 'mo', 'HitTest', 'off');
        gui.pre_cdf_amp_threshold = plot(-1*[1 1], [0 1], 'k', 'HitTest', 'off');
        hold off;
        xlim([0 1]);
        ylim([0 1]);
        xticks(0:0.1:1);
        yticks(0:0.1:1);
        grid on;
        xlabel('Pre-amp / max(pre-amps)');
        ylabel('Pre-amp CDF');
        
        % Setup the LOCAL trace plot
        %------------------------------------------------------------
        if ~isempty(M) % Show movie
            gui.movie = subplot(2,4,8);
            gui.movie_frame = imagesc(ds.cells(cell_idx).im,...
                movie_clim);
            axis image;
            colormap gray;
            hold on;
            boundary = ds.cells(cell_idx).boundary;
            plot(boundary(:,1), boundary(:,2), 'c', 'LineWidth', 2, 'HitTest', 'off');
            com = ds.cells(cell_idx).com;
            plot(com(1), com(2), 'b.');
            
            % Draw nearby cells with fluorescence correlation values.
            % Similar to behavior in 'classify_cells'
            num_neighbors = min(10, ds.num_classified_cells-1);
            gui.neighbor_handles = zeros(num_neighbors, 2); % [Boundary Text]
            
            corr_colors = redblue(201);
            neighbor_inds = ds.get_nearest_sources(cell_idx, num_neighbors);
            for k = 1:num_neighbors
                n = neighbor_inds(k);
                nboundary = ds.cells(n).boundary;
                ncom = ds.cells(n).com;
                ncorr = ds.trace_corrs(cell_idx, n);
                
                % By default, show a thin white boundary of nearby cells
                plot(nboundary(:,1), nboundary(:,2), 'w');
                
                % "Correlation" overlay -- can be toggled
                neighbor_desc = sprintf('%d\n(%.4f)', n, ncorr);
                ncorr = round(ncorr*100) + 101;
                color = corr_colors(ncorr,:); % Convert color to redblue
                gui.neighbor_handles(k,1) = plot(nboundary(:,1), nboundary(:,2),...
                     'Color', color,...
                     'LineWidth',2);
                
                gui.neighbor_handles(k,2) = text(ncom(1), ncom(2), neighbor_desc,...
                     'Color', 'w',...
                     'HorizontalAlignment', 'center',...
                     'FontWeight', 'bold', 'Clipping', 'on');
            end
            hold off;
            
            [height, width, ~] = size(M);
            zoom_half_width = min([height, width])/15;
            xlim(com(1) + zoom_half_width*[-1 1]);
            ylim(com(2) + zoom_half_width*[-1 1]);
            
            gui.local = subplot(2,4,5:7);
        else
            gui.local = subplot(2,1,2);
            
            gui.movie = [];
            gui.movie_frame = [];
            gui.neighbor_handles = [];
        end
        gui.local_orig = plot(trace_orig, 'Color', 0.6*[1 1 1], 'HitTest', 'off');
        hold on;
        gui.local_trace = plot(trace, 'k', 'HitTest', 'off');
        gui.local_dots = plot(trace, 'k.', 'HitTest', 'off');
        trial_starts = ds.trial_indices(:,1);
        gui.local_trials = plot(trial_starts, trace(trial_starts), 'ko', 'HitTest', 'off');
        text_y = trace_display_range(1) + 0.05*diff(trace_display_range);
        gui.local_trials_text = cell(ds.num_trials,1);
        for k = 1:ds.num_trials
            trial_start_k = double(trial_starts(k));
            gui.local_trials_text{k} = text(trial_start_k, text_y,...
                 sprintf('Trial %d', k),...
                 'HorizontalAlignment', 'center',...
                 'HitTest', 'off', 'Clipping', 'on');
        end
        gui.local_cursor_dot = plot(-Inf,trace(1),'ro',...
            'MarkerFaceColor','r',...
            'MarkerSize',6,'HitTest','off');
        gui.local_cursor_bar = plot(-Inf*[1 1], trace_display_range, 'k--', 'HitTest', 'off');
        gui.local_thresh = plot([1 num_frames], -Inf*[1 1], 'm--', 'HitTest', 'off');
        gui.local_auto = plot(-1, -1, 'm:', 'HitTest', 'off');
        gui.local_sel_event = rectangle('Position',[-1 trace_display_range(1) 0 diff(trace_display_range)],...
                  'EdgeColor', 'none',...
                  'FaceColor', [1 0 1 0.2],... // Transparent magenta
                  'HitTest', 'off');
        gui.local_sel_event_peak = plot([-1 -1], trace_display_range, 'm', 'HitTest', 'off');
        gui.local_auto_amps = plot(-1, -1, 'm', 'LineWidth', 2, 'HitTest', 'off');
        gui.local_manual = plot(-1, -1, 'r');
        hold off;
        ylim(trace_display_range);
        xlabel('Frame');
        ylabel('Fluorescence');
        
        % Render data
        %------------------------------------------------------------
        draw_histogram(gui, stats);
        draw_thresholds(gui);
        
        % Add GUI event listeners
        set(gui.global, 'ButtonDownFcn', {@global_plot_handler, gui});
        set(gui.histogram, 'ButtonDownFcn', {@histogram_handler, gui});
        set(gui.post_amps, 'ButtonDownFcn', {@post_cdf_handler, gui});
        set(gui.pre_amps, 'ButtonDownFcn', {@pre_cdf_handler, gui});
        set(gui.local, 'ButtonDownFcn', {@local_plot_handler, gui});
        set(hfig, 'WindowButtonMotionFcn', {@track_cursor, gui});
        set(hfig, 'WindowScrollWheelFcn', {@scroll_plot, gui});
        
        function track_cursor(~, e, gui)
            x = round(e.IntersectionPoint(1));
            if ((1<=x)&&(x<=gui.num_frames))
                if ~isempty(M)
                    set(gui.movie_frame, 'CData', M(:,:,x));
                end
                set(gui.local_cursor_bar,'XData',x*[1 1]);
                set(gui.local_cursor_dot,'XData',x,'YData',trace(x));
            end
        end % track_cursor
        
        function scroll_plot(~, e, gui)
            if (state.show_trials && (ds.num_trials > 1)) % Scroll by trials
                trial_idx = state.last_requested_trial;
                if (e.VerticalScrollCount < 0) % Scroll up
                    trial_idx = trial_idx - 1;
                else
                    trial_idx = trial_idx + 1;
                end
                
                trial_idx = max(1, trial_idx); % Clamp
                trial_idx = min(trial_idx, ds.num_trials);
                show_trial(trial_idx, gui);
                
            else % Default scrolling
                if (e.VerticalScrollCount < 0) % Scroll up
                    get_prev_page(gui);
                else
                    get_next_page(gui);
                end
            end
        end % scroll_plot
        
    end % setup_gui

    % GUI mouse event handlers
    %------------------------------------------------------------
    function global_plot_handler(~, e, gui)
        switch e.Button
            case 1 % Left click -- Move the local viewpoint
                x = round(e.IntersectionPoint(1));
                if ((1<=x) && (x<=gui.num_frames))
                    if state.show_trials
                        t = find(x >= ds.trial_indices(:,1), 1, 'last');
                        show_trial(t, gui);
                    else
                        state.x_anchor = x - state.x_range/2;
                        draw_local_window(gui, state);
                    end
                else
                    fprintf('\n  Not a valid frame for this trace!\n');
                end
                
            case 3 % Right click -- Set threshold
                t = e.IntersectionPoint(2);
                set_thresholds(t, [], [], gui);
        end
    end % global_plot_handler

    function histogram_handler(~, e, gui)
        switch e.Button
            case 1 % Left click
                
            case 3 % Right click -- Set threshold
                t = e.IntersectionPoint(1);
                set_thresholds(t, [], [], gui);
        end
    end % histogram_handler

    function post_cdf_handler(~, e, gui)
        switch e.Button
            case 1 % Left click
                
            case 3 % Right click -- Set the merge threshold
                x = e.IntersectionPoint(1);
                set_thresholds([], x, [], gui);
        end
    end % post_cdf_handler
    
    function pre_cdf_handler(~, e, gui)
        switch e.Button
            case 1 % Left click -- Select a particular event
                x = e.IntersectionPoint(1);
                
                % Find the event with the nearest PRE amplitude
                if ~isempty(events.data)
                    pre_event_amps = events.data(:,3);
                    delta_amp = abs(pre_event_amps - max(pre_event_amps)*x);
                    [~, se] = min(delta_amp);

                    select_event(se, gui);

                    % Move the local window
                    sel_frame = events.data(se,2);
                    state.x_anchor = sel_frame - 1/2 * state.x_range;             
                    draw_local_window(gui, state);
                end
            case 3 % Right click -- Set the amplitude threshold
                x = e.IntersectionPoint(1);
                set_thresholds([], [], x, gui);
        end
    end % pre_cdf_handler

    function local_plot_handler(~, e, gui)
        switch e.Button
            case 1 % Left click
                x = round(e.IntersectionPoint(1));

                % Find the nearest event
                if ~isempty(events.data)
                    event_times = events.data(:,2);
                    delta_times = abs(event_times - x);
                    [~, se] = min(delta_times);

                    select_event(se, gui);
                end
                    
            case 3 % Right click
                
        end
    end % local_plot_handler

    % GUI
    %------------------------------------------------------------   
    function update_gui_state(gui, state)
        draw_local_window(gui, state);
        
        if state.show_orig
            set(gui.local_orig, 'Visible', 'on');
        else
            set(gui.local_orig, 'Visible', 'off');
        end
        
        if state.show_dots
            set(gui.local_dots, 'Visible', 'on');
        else
            set(gui.local_dots, 'Visible', 'off');
        end
        
        if state.show_trials
            vis_val = 'on';
        else
            vis_val = 'off';
        end
        set(gui.local_trials, 'Visible', vis_val);
        for k = 1:ds.num_trials
            set(gui.local_trials_text{k}, 'Visible', vis_val);
        end
        
        if state.show_neighbors
            vis_val = 'on';
        else
            vis_val = 'off';
        end
        num_neighbors = size(gui.neighbor_handles, 1);
        for k = 1:num_neighbors
            set(gui.neighbor_handles(k,1), 'Visible', vis_val);
            set(gui.neighbor_handles(k,2), 'Visible', vis_val);
        end
    end

    function get_next_page(gui)
        current_end = state.x_anchor + state.x_range;
        if (current_end >= gui.num_frames)
            new_anchor = gui.num_frames - state.x_range + 1;
        else
            new_anchor = state.x_anchor + 0.1*state.x_range + 1;
        end

        state.x_anchor = new_anchor;
        draw_local_window(gui, state);
    end % get_next_page

    function get_prev_page(gui)
        new_anchor = state.x_anchor - (0.1*state.x_range + 1);
        state.x_anchor = max(1, new_anchor);
        draw_local_window(gui, state);
    end % get_prev_page

    function draw_local_window(gui, state)
        figure(gui.hfig);
        
        rect_pos = get(gui.global_rect, 'Position');
        rect_pos(1) = state.x_anchor;
        rect_pos(3) = state.x_range;
        set(gui.global_rect, 'Position', rect_pos);
               
        subplot(gui.local);
        xlim([state.x_anchor, state.x_anchor+state.x_range-1]);
    end

    function update_traces(gui)
        set(gui.global_trace, 'YData', trace);
        set(gui.local_trace, 'YData', trace);
        set(gui.local_dots, 'YData', trace);
        set(gui.local_trials, 'YData', trace);
    end

    function draw_histogram(gui, stats)
        subplot(gui.histogram);
        
        set(gui.histogram_data, 'XData', stats.hist_centers,...
                                'YData', stats.hist_counts);

        % First power of 10 that exceeds the maximum count
        count_range = [1 10^ceil(log10(max(stats.hist_counts)))];
        ylim(count_range);
        
        set(gui.histogram_mode, 'XData', events.info.baseline*[1 1],...
                                'YData', count_range);
        set(gui.histogram_thresh, 'XData', events.info.threshold*[1 1],...
                                  'YData', count_range);
    end
    
    function draw_thresholds(gui)
        % Clear previously selected events
        deselect_event(gui);
        
        if ~isempty(events.data)
            event_peak_times = events.data(:,2);
            num_events = length(event_peak_times);
            
            post_amps = sort(events.data(:,4));
            post_cdf_x = post_amps / (max(trace) - events.info.baseline);
            post_cdf_y = (1:num_events)/num_events;
            
            event_amps = sort(events.data(:,3));
            pre_cdf_x = event_amps / event_amps(end);
            pre_cdf_y = (1:num_events)/num_events;
        else
            event_peak_times = [];
            num_events = 0;
            post_cdf_x = [];
            post_cdf_y = [];
            pre_cdf_x = [];
            pre_cdf_y = [];
        end
        
        % GLOBAL subplot
        set(gui.global_thresh, 'YData', events.info.threshold*[1 1]);
        set(gui.global_auto, 'XData', event_peak_times, 'YData', trace(event_peak_times));
        update_global_title(gui);
        
        % HISTOGRAM subplot
        set(gui.histogram_thresh, 'XData', events.info.threshold*[1 1]);
        title(gui.histogram, sprintf('threshold=%.1f', events.info.threshold));
        
        % POST-trough amplitudes CDF subplot
        set(gui.post_cdf_amp_threshold, 'XData', events.info.merge_threshold*[1 1]);
        set(gui.post_cdf, 'XData', post_cdf_x, 'YData', post_cdf_y);
        title(gui.post_amps, sprintf('merge\\_threshold=%.3f', events.info.merge_threshold));
                
        % PRE-trough amplitudes CDF subplot
        set(gui.pre_cdf_amp_threshold, 'XData', events.info.amp_threshold*[1 1]);
        set(gui.pre_cdf, 'XData', pre_cdf_x, 'YData', pre_cdf_y);
        title(gui.pre_amps, sprintf('amp\\_threshold=%.3f', events.info.amp_threshold));
        
        % LOCAL subplot
        set(gui.local_thresh, 'YData', events.info.threshold*[1 1]);
        
        % Note: NaN's break connections between line segments
        X = kron(event_peak_times', [1 1 NaN]);
        Y = repmat([gui.trace_display_range NaN], 1, num_events);
        set(gui.local_auto, 'XData', X, 'YData', Y);
        
        % Draw event amplitudes
        if (num_events > 0)
            Y = zeros(3, num_events);
            Y(1,:) = trace(event_peak_times);
            Y(2,:) = Y(1,:) - events.data(:,3)'; % Peak minus amplitude
            Y(3,:) = NaN;
            Y = Y(:);
        else
            Y = [];
        end
        set(gui.local_auto_amps, 'XData', X, 'YData', Y);
    end % draw_thresholds

    function update_global_title(gui)
        num_events = size(events.data,1);
        
        subplot(gui.global);
        if num_events == 1
            event_str = 'event';
        else
            event_str = 'events';
        end
        title(sprintf('Cell %d, Cutoff freq = %.1f Hz\n%d %s',...
            cell_idx, events.info.cutoff_freq, num_events, event_str));
    end % update_event_tally

    function show_trial(trial_idx, gui)
        trial_start = ds.trial_indices(trial_idx,1);
        trial_end = ds.trial_indices(trial_idx,4);
        trial_range = trial_end - trial_start + 1;

        state.last_requested_trial = trial_idx;
        state.x_anchor = trial_start - 0.25 * trial_range;
        state.x_range = 1.5 * trial_range;
        draw_local_window(gui, state);
    end % show_trial

    function select_event(event_idx, gui)
        % Event index refers to the row of 'events.data'.
        % An index of "0" corresponds to not selecting any event.
        num_events = size(events.data, 1);
        if (0 <= event_idx) && (event_idx <= num_events)
            state.sel_event = event_idx;
            if (event_idx > 0)
                sel_event_onset_frame = events.data(event_idx,1);
                sel_event_peak_frame = events.data(event_idx,2);
                sel_event_duration = sel_event_peak_frame - sel_event_onset_frame;
                
                event_pre_amps = events.data(:,3);
                event_pre_amp = event_pre_amps(event_idx);
                               
                sel_event_pre_normamp = event_pre_amp/max(event_pre_amps);
                sel_event_pre_cdf = sum(event_pre_amp>=event_pre_amps) / num_events;
            else % event_idx == 0
                sel_event_onset_frame = -1;
                sel_event_peak_frame = -Inf;
                sel_event_duration = 0;

                sel_event_pre_normamp = -Inf;
                sel_event_pre_cdf = -1;
            end
            set(gui.pre_cdf_sel_event, 'XData', sel_event_pre_normamp, 'YData', sel_event_pre_cdf);
            rect_pos = get(gui.global_rect, 'Position');
            rect_pos(1) = sel_event_onset_frame;
            rect_pos(3) = sel_event_duration;
            set(gui.local_sel_event, 'Position', rect_pos);
            set(gui.local_sel_event_peak', 'Xdata', sel_event_peak_frame*[1 1]);
        end
    end % select_event

    function deselect_event(gui)
        state.sel_event = 0;
        select_event(0, gui);
    end

    % Data processing
    %------------------------------------------------------------
    function [tr, stats] = preprocess_trace(fps, cutoff_freq)
        tr = filter_trace(ds, cell_idx, fps, cutoff_freq);
        [baseline, sigma, stats] = estimate_baseline_sigma(tr);
        
        events.info.cutoff_freq = cutoff_freq;
        events.info.baseline = baseline;
        events.info.sigma = sigma;
    end

    function compute_events()
        % Core event computation. Uses the parameters in events.info, and
        % fills out events.data
        events.data = find_events_in_trials(trace, ds.trial_indices,...
            events.info.threshold,...
            events.info.baseline,...
            events.info.merge_threshold,...
            events.info.amp_threshold);
    end

    function set_thresholds(threshold, merge_threshold, amp_threshold, gui)
        if ~isempty(threshold)
            events.info.threshold = threshold;
        end
        if ~isempty(merge_threshold)
            merge_threshold = max(0, merge_threshold);
            merge_threshold = min(merge_threshold, 1);
            events.info.merge_threshold = merge_threshold;
        end
        if ~isempty(amp_threshold)
            amp_threshold = max(0, amp_threshold);
            amp_threshold = min(amp_threshold, 1);
            events.info.amp_threshold = amp_threshold;
        end
        compute_events();
        draw_thresholds(gui);
    end % set_threshold

end % detect_events_interactively

function display_range = compute_display_range(trace)
    display_range = [min(trace) max(trace)];
    display_range = display_range + 0.1*diff(display_range)*[-1 1];
end % compute_display_range

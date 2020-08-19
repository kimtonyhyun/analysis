function events = detect_events(ds, cell_idx, varargin)
% Detect calcium events in the fluorescence trace interactively. Core
% computations are:
%   - Determine default thresholds (via 'estimate_baseline_sigma.m'),
%   - Compute events (via 'find_events_in_trials.m').
%
% The application allows for real-time visualization of detected events as
% a function of the algorithmic parameters.

use_prompt = true;
M = [];
fps = 30;
cutoff_freq = 4; % Hz. From Wagner et al. 2019

for j = 1:length(varargin)
    vararg = varargin{j};
    if ischar(vararg)
        switch lower(vararg)
            case 'noprompt'
                use_prompt = false;
            case 'fps'
                fps = varargin{j+1};
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

% Default event processing
%------------------------------------------------------------

% We look for events in a smoothed version of the trace
trace_orig = ds.get_trace(cell_idx);
trace = filter_trace(ds, cell_idx, fps, cutoff_freq);
num_frames = length(trace);

% Default threshold values
[baseline, sigma, stats] = estimate_baseline_sigma(trace);

init_info = struct('num_frames', num_frames,...
                   'baseline', baseline,...
                   'sigma', sigma,...
                   'threshold', baseline + 5*sigma,...
                   'amp_threshold', 0.1);
events = struct('info', init_info, 'data', []);
compute_events(); % Fills in 'events.data'

% Set up GUI
%------------------------------------------------------------
if ~use_prompt
    return;
end

if ~isempty(M) && isempty(movie_clim)
    movie_clim = compute_movie_scale(M);
end

hfig = figure;
gui = setup_gui(hfig, num_frames, compute_display_range(trace), stats, trace_orig);

% GUI state
state.x_anchor = 1;
state.x_range = min(500, num_frames);
state.show_orig = true;
state.show_dots = false;
state.show_trials = (ds.num_trials > 1);
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
    resp = lower(strtrim(input(prompt, 's')));
    val = str2double(resp);

    if (~isnan(val)) % Is a number
        if (1 <= val) && (val <= ds.num_trials)
            show_trial(val, gui);
        end
    else % Not a number
        switch (resp)
            case 'q' % "quit"
                close(hfig);
                break;

            case {'i', 'z'} % zoom in
                x_center = state.x_anchor + 1/2 * state.x_range;
                state.x_range = 0.5*state.x_range;
                state.x_anchor = x_center - 1/2 * state.x_range;
                redraw_local_window(gui, state);
                
            case 'o' % zoom out
                x_center = state.x_anchor + 1/2 * state.x_range;
                state.x_range = 2*state.x_range;
                state.x_anchor = x_center - 1/2 * state.x_range;
                redraw_local_window(gui, state);
                
            case 'r' % toggle "raw"/original trace
                state.show_orig = ~state.show_orig;
                update_gui_state(gui, state);
                
            case 'd' % show dots
                state.show_dots = ~state.show_dots;
                update_gui_state(gui, state);
                
            case 'b' % show trials
                state.show_trials = ~state.show_trials;
                update_gui_state(gui, state);
                
            case {'', 'n'} % next trial
                if state.last_requested_trial < ds.num_trials
                    show_trial(state.last_requested_trial+1, gui);
                end
                
            case 't' % reset threshold
                events.info = init_info;
                compute_events();
                draw_thresholds(gui);

            otherwise
                fprintf('  Sorry, could not parse "%s"\n', resp);
        end
    end
end % Main interaction loop

    % Supplementary functions
    %------------------------------------------------------------
    function gui = setup_gui(hf, num_frames, trace_display_range, stats, trace_orig)
        % Display parameters kept around for convenience
        gui.hfig = hf;
        gui.num_frames = num_frames;
        gui.trace_display_range = trace_display_range;
        
        % Setup the GLOBAL trace plot
        %------------------------------------------------------------
        gui.global = subplot(2,4,1:3);
        gui.global_rect = rectangle('Position',[-1 trace_display_range(1) 1 diff(trace_display_range)],...
                  'EdgeColor', 'none',...
                  'FaceColor', 'c', 'HitTest', 'off');
        hold on;
        plot(trace, 'k', 'HitTest', 'off');
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
        gui.histogram = subplot(4,4,4);
        semilogy(stats.hist_centers, stats.hist_counts, 'k.', 'HitTest', 'off');
        xlim(trace_display_range);
        hold on;
        % First power of 10 that exceeds the maximum count
        count_range = [1 10^ceil(log10(max(stats.hist_counts)))];
        ylim(count_range);
        for k = 1:size(stats.percentiles,1)
            f = stats.percentiles(k,2);
            plot(f*[1 1], count_range, 'Color', 0.5*[1 1 1], 'HitTest', 'off');
        end
        gui.histogram_thresh = plot(-Inf*[1 1], count_range, 'm', 'HitTest', 'off');
        hold off;
        xlabel('Fluorescence');
        ylabel('Trace histogram');

        % Setup the event amplitude CDF
        %------------------------------------------------------------
        gui.event_amp_cdf = subplot(4,4,8);
        gui.cdf = plot(-1, -1, 'm.-', 'HitTest', 'off');
        hold on;
        gui.cdf_sel_event = plot(-1, -1, 'mo', 'HitTest', 'off');
        gui.cdf_amp_threshold = plot(-1*[1 1], [0 1], 'k', 'HitTest', 'off');
        hold off;
        xlim([0 1]);
        ylim([0 1]);
        xticks(0:0.1:1);
        yticks(0:0.1:1);
        grid on;
        xlabel('Norm event amplitude');
        ylabel('CDF');
        
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
            num_neighbors = min(10, ds.num_classified_cells-1);
            neighbor_inds = ds.get_nearest_sources(cell_idx, num_neighbors);
            for n = neighbor_inds
                boundary = ds.cells(n).boundary;
                plot(boundary(:,1), boundary(:,2), 'w');
                
                ncom = ds.cells(n).com;
                text(ncom(1), ncom(2), num2str(n),...
                     'Color', 'w',...
                     'HorizontalAlignment', 'center',...
                     'HitTest', 'off', 'Clipping', 'on');
            end
            hold off;
            
            [height, width, ~] = size(M);
            zoom_half_width = min([height, width])/15;
            xlim(com(1) + zoom_half_width*[-1 1]);
            ylim(com(2) + zoom_half_width*[-1 1]);
            
            gui.local = subplot(2,4,5:7);
        else
            gui.local = subplot(2,1,2);
        end
        gui.local_orig = plot(trace_orig, 'Color', 0.6*[1 1 1], 'HitTest', 'off');
        hold on;
        plot(trace, 'k', 'HitTest', 'off');
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
        
        draw_thresholds(gui);
        
        % Add GUI event listeners
        set(gui.global, 'ButtonDownFcn', {@global_plot_handler, gui});
        set(gui.histogram, 'ButtonDownFcn', {@histogram_handler, gui});
        set(gui.event_amp_cdf, 'ButtonDownFcn', {@cdf_handler, gui});
        set(gui.local, 'ButtonDownFcn', {@local_plot_handler, gui});
        set(hf, 'WindowButtonMotionFcn', {@track_cursor, gui});
        set(hf, 'WindowScrollWheelFcn', {@scroll_plot, gui});
        
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
                        redraw_local_window(gui, state);
                    end
                else
                    fprintf('\n  Not a valid frame for this trace!\n');
                end
                
            case 3 % Right click -- Set threshold
                t = e.IntersectionPoint(2);
                set_thresholds(t, [], gui);
        end
    end % global_plot_handler

    function histogram_handler(~, e, gui)
        switch e.Button
            case 1 % Left click
                
            case 3 % Right click -- Set threshold
                t = e.IntersectionPoint(1);
                set_thresholds(t, [], gui);
        end
    end % histogram_handler

    function cdf_handler(~, e, gui)
        switch e.Button
            case 1 % Left click -- Select a particular event
                x = e.IntersectionPoint(1);
                
                % Find the event with the nearest amplitude
                event_amps = events.data(:,3);
                delta_amp = abs(event_amps - max(event_amps)*x);
                [~, se] = min(delta_amp);
                
                show_event(se, gui);
                
                % Move the local window
                sel_frame = events.data(se,2);
                state.x_anchor = sel_frame - 1/2 * state.x_range;             
                redraw_local_window(gui, state);
            case 3 % Right click -- Set the amplitude threshold
                x = e.IntersectionPoint(1);
                set_thresholds([], x, gui);
        end
    end % cdf_handler

    function local_plot_handler(~, e, gui)
        switch e.Button
            case 1 % Left click
                x = round(e.IntersectionPoint(1));

                % Find the nearest event
                event_times = events.data(:,2);
                delta_times = abs(event_times - x);
                [~, se] = min(delta_times);

                show_event(se, gui);
                    
            case 3 % Right click
                
        end
    end % local_plot_handler

    % GUI
    %------------------------------------------------------------   
    function update_gui_state(gui, state)
        redraw_local_window(gui, state);
        
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
            trials_vis = 'on';
        else
            trials_vis = 'off';
        end
        set(gui.local_trials, 'Visible', trials_vis);
        for k = 1:ds.num_trials
            set(gui.local_trials_text{k}, 'Visible', trials_vis);
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
        redraw_local_window(gui, state);
    end % get_next_page

    function get_prev_page(gui)
        new_anchor = state.x_anchor - (0.1*state.x_range + 1);
        state.x_anchor = max(1, new_anchor);
        redraw_local_window(gui, state);
    end % get_prev_page

    function redraw_local_window(gui, state)
        figure(gui.hfig);
        
        rect_pos = get(gui.global_rect, 'Position');
        rect_pos(1) = state.x_anchor;
        rect_pos(3) = state.x_range;
        set(gui.global_rect, 'Position', rect_pos);
               
        subplot(gui.local);
        xlim([state.x_anchor, state.x_anchor+state.x_range-1]);
    end
    
    function draw_thresholds(gui)
        % Clear previously selected events
        clear_event(gui);
        
        if ~isempty(events.data)
            event_peak_times = events.data(:,2);
            num_events = length(event_peak_times);
            
            % Assumes events are sorted by amplitude
            event_amps = sort(events.data(:,3));
            cdf_x = event_amps / event_amps(end);
            cdf_y = (1:num_events)/num_events;
        else
            event_peak_times = [];
            num_events = 0;
            cdf_x = [];
            cdf_y = [];
        end
        
        % GLOBAL subplot
        set(gui.global_thresh, 'YData', events.info.threshold*[1 1]);
        set(gui.global_auto, 'XData', event_peak_times, 'YData', trace(event_peak_times));
        update_event_tally(gui);
        
        % HISTOGRAM subplot
        set(gui.histogram_thresh, 'XData', events.info.threshold*[1 1]);
        
        % CDF subplot
        set(gui.cdf_amp_threshold, 'XData', events.info.amp_threshold*[1 1]);
        set(gui.cdf, 'XData', cdf_x, 'YData', cdf_y);
        
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

    function update_event_tally(gui)
        num_events = size(events.data,1);
        
        subplot(gui.global);
        title(sprintf('Cell %d: %d events', cell_idx, num_events));
    end % update_event_tally

    function show_trial(trial_idx, gui)
        trial_start = ds.trial_indices(trial_idx,1);
        trial_end = ds.trial_indices(trial_idx,4);
        trial_range = trial_end - trial_start + 1;

        state.last_requested_trial = trial_idx;
        state.x_anchor = trial_start - 0.25 * trial_range;
        state.x_range = 1.5 * trial_range;
        redraw_local_window(gui, state);
    end % show_trial

    function show_event(event_idx, gui)
        % Event index refers to the row of 'events.data'.
        % An index of "0" corresponds to not selecting any event.
        num_events = size(events.data, 1);
        if (0 <= event_idx) && (event_idx <= num_events)
            state.sel_event = event_idx;
            if (event_idx > 0)
                sel_event_onset_frame = events.data(event_idx,1);
                sel_event_peak_frame = events.data(event_idx,2);
                sel_event_duration = sel_event_peak_frame - sel_event_onset_frame;
                
                event_amps = events.data(:,3);
                event_amp = event_amps(event_idx);
                               
                sel_event_normamp = event_amp/max(event_amps);
                sel_event_cdf = sum(event_amp>=event_amps) / num_events;
            else % event_idx == 0
                sel_event_onset_frame = -1;
                sel_event_peak_frame = -Inf;
                sel_event_duration = 0;
                sel_event_normamp = -Inf;
                sel_event_cdf = -1;
            end
            set(gui.cdf_sel_event, 'XData', sel_event_normamp, 'YData', sel_event_cdf);
            rect_pos = get(gui.global_rect, 'Position');
            rect_pos(1) = sel_event_onset_frame;
            rect_pos(3) = sel_event_duration;
            set(gui.local_sel_event, 'Position', rect_pos);
            set(gui.local_sel_event_peak', 'Xdata', sel_event_peak_frame*[1 1]);
        end
    end % show_event

    function clear_event(gui)
        show_event(0, gui);
    end

    % Data processing
    %------------------------------------------------------------
    function compute_events()
        % Core event computation. Uses the parameters in events.info, and
        % fills out events.data
        events.data = find_events_in_trials(trace, ds.trial_indices,...
            events.info.threshold,...
            events.info.baseline,...
            events.info.amp_threshold);
    end

    function set_thresholds(threshold, amp_threshold, gui)
        if ~isempty(threshold)
            events.info.threshold = threshold;
        end
        if ~isempty(amp_threshold)
            % Clamp to [0 1]
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

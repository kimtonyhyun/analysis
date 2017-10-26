function events = detect_events_interactively(ds, cell_idx, varargin)

use_filter = true;
M = [];
movie_clim = [];
fps = 30;
cutoff_freq = [];

for i = 1:length(varargin)
    vararg = varargin{i};
    if ischar(vararg)
        switch lower(vararg)
            case 'fps'
                fps = varargin{i+1};
            case 'cutoff'
                cutoff_freq = varargin{i+1};
            case 'nofilter'
                use_filter = false;
            case {'m', 'movie'}
                M = varargin{i+1};
            case {'clim', 'movie_clim'}
                movie_clim = varargin{i+1};
        end
    end
end

trace_orig = ds.get_trace(cell_idx);
if ~isempty(M) && isempty(movie_clim)
    movie_clim = compute_movie_scale(M);
end

% We'll for events in a smoothed version of the trace
% Default parameters comes from cerebellar processing, where we used
%   - 30 Hz sampling frequency
%   - 4 Hz cutoff frequency
if isempty(cutoff_freq)
    cutoff_freq = 4/30 * fps;
end
if use_filter
    fprintf('Applying LPF (fc=%.1f Hz) to trace...\n', cutoff_freq);
%     trace = filter_trace(trace_orig, cutoff_freq, fps);

    filt_traces = cell(1,ds.num_trials);
    for tidx = 1:ds.num_trials
        filt_traces{tidx} = filter_trace(ds.get_trace(cell_idx,tidx), cutoff_freq, fps);
    end
    trace = cell2mat(filt_traces);
else
    trace = trace_orig;
end

% Basic trace properties
num_frames = length(trace);
trace_display_range = [min(trace) max(trace)];
trace_display_range = trace_display_range + 0.1*diff(trace_display_range)*[-1 1];
[init_threshold, stats] = estimate_baseline_threshold(trace);

% Application state
state.allow_manual_events = false;
state.x_anchor = 1;
state.x_range = min(500, num_frames);
state.show_orig = true;
state.show_dots = false;
state.show_trials = (ds.num_trials > 1);
state.sel_event = 0;
state.last_requested_trial = 0;

events = struct('threshold', [], 'auto', [], 'manual', []);

hfig = figure;
gui = setup_gui(hfig, num_frames, trace_display_range, stats, trace_orig);
set_threshold(init_threshold, gui);
update_gui_state(gui, state);

% Interaction loop:
%------------------------------------------------------------
prompt = '  >> ';
resp = lower(strtrim(input(prompt, 's')));
val = str2double(resp);

while (1)
    if (~isnan(val)) % Is a number
        if (1 <= val) && (val <= ds.num_trials)
            set_trial(val, gui);
        end
    else % Not a number
        switch (resp)
            case 'q' % "quit"
                close(hfig);
                break;

            case 'z' % zoom in
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
                    set_trial(state.last_requested_trial+1, gui);
                end
                
            case 't' % reset threshold
                set_threshold(init_threshold, gui);
                                
%             case 'm' % toggle manual input -- DISABLED for now
%                 state.allow_manual_events = ~state.allow_manual_events;
                
%             case 'x' % erase last event
%                 num_events = length(events.manual);
%                 if (num_events > 0)
%                     events.manual = events.manual(1:num_events-1);
%                     redraw_manual_events(gui);
%                 else
%                     fprintf('  No events to remove!\n');
%                 end

            otherwise
                fprintf('  Sorry, could not parse "%s"\n', resp);
        end
    end
        
    resp = lower(strtrim(input(prompt, 's')));
    val = str2double(resp);
end % Main interaction loop

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

    function gui = setup_gui(hf, num_frames, trace_display_range, stats, trace_orig)
        % Display parameters kept around for convenience
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
        gui.local_sel_event = plot([-1 -1], trace_display_range, 'm', 'HitTest', 'off');
        gui.local_auto_amps = plot(-1, -1, 'm', 'LineWidth', 2, 'HitTest', 'off');
        gui.local_manual = plot(-1, -1, 'r');
        hold off;
        ylim(trace_display_range);
        grid on;
        xlabel('Frame');
        ylabel('Fluorescence');
        
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
                if state.allow_manual_events
                    x = seek_localmax(trace, x);
                end
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
                set_trial(trial_idx, gui);
                
            else % Default scrolling
                if (e.VerticalScrollCount < 0) % Scroll up
                    get_prev_page(gui);
                else
                    get_next_page(gui);
                end
            end
        end % scroll_plot
        
    end % setup_gui

    % Update the GUI
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

    function redraw_local_window(gui, state)
        rect_pos = get(gui.global_rect, 'Position');
        rect_pos(1) = state.x_anchor;
        rect_pos(3) = state.x_range;
        set(gui.global_rect, 'Position', rect_pos);
               
        subplot(gui.local);
        xlim([state.x_anchor, state.x_anchor+state.x_range-1]);
    end
    
    function redraw_threshold(gui)
        % GLOBAL subplot
        set(gui.global_thresh, 'YData', events.threshold*[1 1]);
        auto_peaks = events.auto(:,2)';
        set(gui.global_auto, 'XData', auto_peaks, 'YData', trace(auto_peaks));
        update_event_tally(gui);
        
        % HISTOGRAM subplot
        set(gui.histogram_thresh, 'XData', events.threshold*[1 1]);
        
        % CDF subplot
        num_auto_events = size(events.auto, 1);
        max_amp = events.auto(end,3);
        set(gui.cdf, 'XData', events.auto(:,3)/max_amp,...
                     'YData', (1:num_auto_events)/num_auto_events);
        
        % LOCAL subplot
        set(gui.local_thresh, 'YData', events.threshold*[1 1]);
        
        % Note: NaN's break connections between line segments
        X = kron(auto_peaks, [1 1 NaN]);
        Y = repmat([gui.trace_display_range NaN], 1, num_auto_events);
        set(gui.local_auto, 'XData', X, 'YData', Y);
        
        % Draw event amplitudes
        Y = zeros(3, num_auto_events);
        Y(1,:) = trace(events.auto(:,2)); % Peak
        Y(2,:) = Y(1,:) - events.auto(:,3)'; % Peak minus amplitude
        Y(3,:) = NaN;
        set(gui.local_auto_amps, 'XData', X, 'YData', Y(:));
    end % redraw_threshold

    function redraw_manual_events(gui)
        set(gui.global_manual, 'XData', events.manual, 'YData', trace(events.manual));
        
        X = kron(events.manual', [1 1 NaN]);
        Y = repmat([gui.trace_display_range NaN], 1, length(events.manual));
        set(gui.local_manual, 'XData', X, 'YData', Y);
        
        update_event_tally(gui);
    end % redraw_manual_events

    function update_event_tally(gui)
        num_auto = size(events.auto,1);
        num_manual = length(events.manual);
        
        subplot(gui.global);
        title(sprintf('Num events: %d (auto), %d (manual)', num_auto, num_manual));
    end % update_event_tally

    % Event handlers for mouse input
    %------------------------------------------------------------
    function global_plot_handler(~, e, gui)
        switch e.Button
            case 1 % Left click -- Move the local viewpoint
                x = round(e.IntersectionPoint(1));
                if ((1<=x) && (x<=gui.num_frames))
                    state.x_anchor = x - state.x_range/2;
                    redraw_local_window(gui, state);
                    
                    % To have consistent behavior with subsequent mouse
                    % scrolls, set last_requested_trial
                    state.last_requested_trial = find(x >= ds.trial_indices(:,1)', 1, 'last');
                else
                    fprintf('\n  Not a valid frame for this trace!\n');
                end
                
            case 3 % Right click -- Set threshold
                t = e.IntersectionPoint(2);
                set_threshold(t, gui);
        end
    end % global_plot_handler

    function histogram_handler(~, e, gui)
        switch e.Button
            case 1 % Left click
                
            case 3 % Right click -- Set threshold
                t = e.IntersectionPoint(1);
                set_threshold(t, gui);
        end
    end % histogram_handler

    function cdf_handler(~, e, gui)
        switch e.Button
            case 1 % Left click -- "Select" a particular event
                x = e.IntersectionPoint(1);
                
                % Find the event with the nearest amplitude
                event_amps = events.auto(:,3);
                delta_amp = abs(event_amps - max(event_amps)*x);
                [~, se] = min(delta_amp);
                
                select_event(se, gui);
                
                % Move the local window
                sel_frame = events.auto(se,2);
                state.x_anchor = sel_frame - 1/2 * state.x_range;             
                redraw_local_window(gui, state);
            case 3 % Right click
                
        end
    end % cdf_handler

    function local_plot_handler(~, e, gui)
        switch e.Button
            case 1 % Left click
                x = round(e.IntersectionPoint(1));
                if state.allow_manual_events
                    add_manual_event(x, gui);
                else
                    % Find the nearest event
                    event_times = events.auto(:,2);
                    delta_times = abs(event_times - x);
                    [~, se] = min(delta_times);
                    
                    select_event(se, gui);
                end
            case 3 % Right click
                
        end
    end % local_plot_handler

    % Data processing
    %------------------------------------------------------------
    function add_manual_event(x, gui)
        if ((1<=x) && (x<=gui.num_frames))
            x = seek_localmax(trace, x);
            % Don't make duplicate events
            auto_peaks = events.auto(:,2);
            if ~ismember(x, auto_peaks) && ~ismember(x, events.manual)
                events.manual = [events.manual; x];
                redraw_manual_events(gui);
            end
        else
            fprintf('\n  Not a valid event for this trace!\n');
        end
    end % add_manual_event

    function set_threshold(t, gui)
        events.threshold = t;
        es = find_events_in_trials(trace, ds.trial_indices, t, stats.mode);
        events.auto = cell2mat(es);
        events.auto = sortrows(events.auto, 3); % Sort events by amplitude
        
        select_event(0, gui);
        redraw_threshold(gui);
    end % set_threshold

    function set_trial(trial_idx, gui)
        trial_start = ds.trial_indices(trial_idx,1);
        trial_end = ds.trial_indices(trial_idx,4);
        trial_range = trial_end - trial_start + 1;

        state.last_requested_trial = trial_idx;
        state.x_anchor = trial_start - 0.25 * trial_range;
        state.x_range = 1.5 * trial_range;
        redraw_local_window(gui, state);
    end % set_trial

    function select_event(event_idx, gui)
        % Event index refers to the row of 'events.auto'. An index of "0"
        % corresponds to not selecting any event.
        num_events = size(events.auto, 1);
        if (0 <= event_idx) && (event_idx <= num_events)
            state.sel_event = event_idx;
            if (event_idx > 0)
                event_amps = events.auto(:,3);
                sel_event_frame = events.auto(event_idx,2);
                sel_event_normamp = event_amps(event_idx)/max(event_amps);
                sel_event_cdf = event_idx / num_events;
            else
                sel_event_frame = -Inf;
                sel_event_normamp = -Inf;
                sel_event_cdf = -1;
            end
            set(gui.cdf_sel_event, 'XData', sel_event_normamp, 'YData', sel_event_cdf);
            set(gui.local_sel_event, 'XData', sel_event_frame*[1 1]);
        end
    end % select_event

end % detect_events_interactively

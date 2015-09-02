function [resp, movie_clim] = view_cell_interactively(ds, cell_idx, movie, fps, movie_clim)
% Visually inspect the active portions of an IC trace side-by-side with
%   the provided miniscope movie.
%
% Do NOT modify the large 'movie' matrix, lest Matlab make a duplicate
%   in memory!
%

filter = ds.cells(cell_idx).im;
trace_orig  = ds.get_trace(cell_idx);
trace_fixed = fix_baseline(trace_orig);
time = 1/fps*((1:length(trace_orig))-1);

% Some parameters
ic_filter_threshold = 0.3; % For generating the IC filter outline
mad_scale = 10; % Used for coarse detection of activity in the IC trace
active_frame_padding = 5*fps; % Use 100 for 20 Hz movie
time_window = 100/fps; % Width of running window

% Generate the outline of the filter
%------------------------------------------------------------
subplot(3,3,[4 5 7 8]);
h = imagesc(rescale_filter_to_clim(filter, movie_clim), movie_clim);
colormap gray;
axis image;
xlabel('x [px]');
ylabel('y [px]');
hold on;
[ic_boundaries, ic_mask] = compute_ic_boundary(filter, ic_filter_threshold);
for i = 1:length(ic_boundaries)
    ic_boundary = ic_boundaries{i};
    plot(ic_boundary(:,1), ic_boundary(:,2), 'c', 'LineWidth', 2);
end

% Plot boundaries of other cells, and retrieve their handles so that
%   we can toggle the boundaries on and off
other_cells = setdiff(1:ds.num_cells, cell_idx);
num_other_cells = length(other_cells);
other_cell_handles = zeros(num_other_cells, 1);

for n = 1:num_other_cells
    oc_idx = other_cells(n);
    boundary = ds.cells(oc_idx).boundary;
    if isempty(ds.cells(oc_idx).label)
        color = 'w';
    else
        if ds.is_cell(oc_idx)
            color = 'g';
        else
            color = 'r';
        end
    end
    other_cell_handles(n) = plot(boundary(:,1), boundary(:,2), color);
end
show_other_cells(false); % Turn off the boundaries of other cells

% Compute the center of mass of the filter
masked_filter = ic_mask .* filter;
[height, width] = size(masked_filter);
COM = [(1:width) * sum(masked_filter,1)';
       (1:height)* sum(masked_filter,2)];
COM = COM / sum(masked_filter(:));
plot(COM(1), COM(2), 'b.');
hold off;

% Start off zoomed
zoom_half_width = min([width, height])/10;
xlim(COM(1)+zoom_half_width*[-1 1]);
ylim(COM(2)+zoom_half_width*[-1 1]);

% Compute the active portions of the trace
%------------------------------------------------------------
trace = trace_orig; % By default, select the original trace
thresh = mad_scale * compute_mad(trace);

[active_periods, num_active_periods] =...
    parse_active_frames(trace > thresh, active_frame_padding);

setup_traces();

% Interaction loop:
%   Display the user-specified active period
%------------------------------------------------------------
prompt = 'Cell viewer >> ';
resp = lower(strtrim(input(prompt, 's')));
val = str2double(resp);

% State of interaction loop
state.last_val = [];
state.zoomed = true;
state.show_other_cells = false;
state.baseline_removed = false;

while (1)
    if (~isnan(val)) % Is a number
        if ((1 <= val) && (val <= num_active_periods))
            display_active_period(val);
            state.last_val = val;
        else
            fprintf('  Sorry, %d is not a valid period index for this IC\n', val);
        end
    else % Not a number
        switch (resp)
            case {'', 'q'} % "quit"
                break;
                
            case 'a' % "all"
                display_active_period(1:num_active_periods);
                
            case 't' % "threshold"
                fprintf('  Please select a new threshold on the global trace\n');
                while (1)
                    [~, thresh] = ginput(1);
                    if (gca == global_trace)
                        break;
                    else
                        fprintf('  Error! New threshold must be defined on the GLOBAL trace\n');
                    end
                end
                fprintf('  New threshold value of %.3f selected!\n', thresh);

                % Recompute active periods and redraw
                [active_periods, num_active_periods] =...
                    parse_active_frames(trace > thresh, active_frame_padding);
                setup_traces();
                
            case 'b' % Fix "baseline"
                if (state.baseline_removed)
                    trace = trace_orig;
                    fprintf('  Showing original trace without baseline correction\n');
                else
                    trace = trace_fixed;
                    fprintf('  Showing trace with baseline correction\n');
                end
                
                thresh = mad_scale * compute_mad(trace);
                [active_periods, num_active_periods] = ...
                    parse_active_frames(trace > thresh, active_frame_padding);
                setup_traces();
                
                state.baseline_removed = ~state.baseline_removed; % Toggle
                
            case 'r' % "replay"
                if ~isempty(state.last_val)
                    display_active_period(state.last_val);
                end
                            
            case 'z' % "zoom"
                subplot(3,3,[4 5 7 8]); % Focus on the movie subplot
                if (state.zoomed) % Return to original view
                    xlim([1 width]);
                    ylim([1 height]);
                    state.zoomed = false;
                else
                    xlim(COM(1)+zoom_half_width*[-1 1]);
                    ylim(COM(2)+zoom_half_width*[-1 1]);
                    state.zoomed = true;
                end
                    
            case {'h', 'l'} % "higher/lower contrast"
                subplot(3,3,[4 5 7 8]); % Focus on the movie subplot
                c_range = diff(movie_clim);
                if (strcmp(resp, 'h'))
                    movie_clim = movie_clim + c_range*[0.1 -0.1];
                    fprintf('  Increased contrast (new CLim=[%.3f %.3f])\n',...
                        movie_clim(1), movie_clim(2));
                else
                    movie_clim = movie_clim + c_range*[-0.1 0.1];
                    fprintf('  Decreased contrast (new CLim=[%.3f %.3f])\n',...
                        movie_clim(1), movie_clim(2));
                end
                set(gca, 'CLim', movie_clim);
                
            case {'n', 'm'} % Show "neighboring" cells
                state.show_other_cells = ~state.show_other_cells;
                show_other_cells(state.show_other_cells);

            otherwise
                fprintf('  Sorry, could not parse "%s"\n', resp);
        end
    end

    resp = lower(strtrim(input(prompt, 's')));
    val = str2double(resp);
end

    % Display subroutines
    %------------------------------------------------------------
    function setup_traces()
        global running_trace t_g t_r dot;
        
        x_range = [time(1) time(end)];
        y_range = [min(trace(:)) max(trace(:))];
        y_delta = y_range(2) - y_range(1);
        y_range = y_range + 0.1*y_delta*[-1 1];
        
        % Prepare global trace
        global_trace = subplot(3,3,[1 2 3]);
        plot(time, trace, 'b');
        hold on;
        plot(x_range, thresh*[1 1], 'r--'); % Display threshold
        for period_idx = 1:num_active_periods
            active_period = active_periods(period_idx, :);
            active_frames = active_period(1):active_period(2);
            plot(time(active_frames), trace(active_frames), 'r');
            text(double(time(active_frames(1))),... % 'text' fails on single
                 double(y_range(2)),...
                 num2str(period_idx),...
                 'Color', 'r',...
                 'VerticalAlignment', 'top');
        end
        xlim(x_range);
        ylim(y_range);
        t_g = plot(time(1)*[1 1], y_range, 'k'); % Time indicator
        xlabel('Time [s]');
        ylabel('Signal [a.u.]');
        hold off;

        % Prepare running trace
        running_trace = subplot(3,3,[6 9]);
        plot(time, trace, 'b');
        hold on;
        for period_idx = 1:num_active_periods
            active_period = active_periods(period_idx, :);
            active_frames = active_period(1):active_period(2);
            plot(time(active_frames), trace(active_frames), 'r');
        end
        xlim([0 time_window]);
        ylim(y_range);
        t_r = plot(time(1)*[1 1], y_range, 'k'); % Time indicator
        dot = plot(time(1), trace(1), 'or',...
                    'MarkerFaceColor', 'r',...
                    'MarkerSize', 12); % Dot
        xlabel('Time [s]');
        ylabel('Signal [a.u.]');
        hold off;
    end % setup_traces
    
    function display_active_period(selected_indices)
        global running_trace t_g t_r dot;
        
        for selected_idx = selected_indices
            frames = active_periods(selected_idx,1):...
                     active_periods(selected_idx,2);
            for k = frames
                A = movie(:,:,k);
                set(h, 'CData', A);

                % Update time indicators and dot
                set(t_g, 'XData', time(k)*[1 1]);
                set(t_r, 'XData', time(k)*[1 1]);
                set(dot, 'XData', time(k), 'YData', trace(k));

                % Update running trace
                set(running_trace, 'XLim', time(k) + time_window/2*[-1 1]);
                drawnow;
            end
        end
    end % display_active_period

    function show_other_cells(show)
        vis_val = 'off';
        if (show)
            vis_val = 'on';
        end
        
        for m = 1:num_other_cells
            set(other_cell_handles(m), 'Visible', vis_val);
        end
    end % show_other_cells

end % main function

function [active_frames, num_active] = parse_active_frames(binary_trace, half_width)
% Segment the active portions of a binary trace into intervals

    if (half_width > 0)
        trace = single(binary_trace);
        trace = conv(trace, ones(1, 2*half_width + 1), 'same');
        trace = logical(trace);
    else
        trace = binary_trace;
    end
    trace_comp = ~trace; % Complement of the trace

    active_frames = [];

    % Loop to find all activity transitions in the trace
    curr = 1;
    while (1)
        next = find(trace(curr:end), 1, 'first');
        if (isempty(next))
            break;
        end
        active_frames = [active_frames curr+(next-1)]; %#ok<*AGROW>
        curr = curr + next;

        next = find(trace_comp(curr:end), 1, 'first');
        if (isempty(next))
            break;
        end
        active_frames = [active_frames curr+(next-2)];
        curr = curr + next;
    end

    if (mod(length(active_frames),2) == 1) % Ended with active frame
        active_frames = [active_frames length(binary_trace)];
    end

    active_frames = reshape(active_frames, 2, length(active_frames)/2)';
    num_active = size(active_frames, 1);
end

function filter_out = rescale_filter_to_clim(filter, clim)
% Numerically rescale the filter to match the provided clim

    f_max = max(filter(:));
    f_min = min(filter(:));

    filter_norm = (filter-f_min)/f_max; % Matched to range [0, 1]

    clim_delta = clim(2)-clim(1);
    clim_usage = 0.8;
    filter_out = clim(1) + (1-clim_usage)/2*clim_delta +...
                 clim_usage*clim_delta*filter_norm;

end
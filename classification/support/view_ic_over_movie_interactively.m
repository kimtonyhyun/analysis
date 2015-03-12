function [resp, movie_clim] = view_ic_over_movie_interactively(ic_filter, time, trace, movie, fps, movie_clim)
% Visually inspect the active portions of an IC trace side-by-side with
%   the provided miniscope movie
%
% 2015 01 31 Tony Hyun Kim

% Some parameters
ic_filter_threshold = 0.3; % For generating the IC filter outline
mad_scale = 8; % Used for coarse detection of activity in the IC trace
active_frame_padding = 5*fps; % Use 100 for 20 Hz movie
time_window = 10; % Width of running window

% Generate the outline of the IC filter
%------------------------------------------------------------
subplot(3,3,[4 5 7 8]);
h = imagesc(ic_filter, movie_clim);
colormap gray;
axis image;
xlabel('x [px]');
ylabel('y [px]');
hold on;
[ic_boundaries, ic_mask] = compute_ic_boundary(ic_filter, ic_filter_threshold);
for i = 1:length(ic_boundaries)
    ic_boundary = ic_boundaries{i};
    plot(ic_boundary(:,1), ic_boundary(:,2), 'r', 'LineWidth', 2);
end

% Compute the center of mass of the filter
masked_filter = ic_mask .* ic_filter;
[height, width] = size(masked_filter);
COM = [(1:width) * sum(masked_filter,1)';
       (1:height)* sum(masked_filter,2)];
COM = COM / sum(masked_filter(:));
plot(COM(1), COM(2), 'b.');
hold off;

zoom_half_width = min([width, height])/10;

% Compute the active portions of the trace
%------------------------------------------------------------
x_range = [time(1) time(end)];
y_range = [min(trace(:)) max(trace(:))];
y_delta = y_range(2) - y_range(1);
y_range = y_range + 0.1*y_delta*[-1 1];

mad = compute_mad(trace);
thresh = mad_scale * mad;

active_periods = parse_active_frames(trace > thresh,...
                                     active_frame_padding);
num_active_periods = size(active_periods, 1);

% Prepare global trace
subplot(3,3,[1 2 3]);
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
t1 = plot(time(1)*[1 1], y_range, 'k'); % Time indicator
xlabel('Time [s]');
ylabel('Signal [a.u.]');
hold off;

% Prepare running trace
a = subplot(3,3,[6 9]);
plot(time, trace, 'b');
hold on;
for period_idx = 1:num_active_periods
    active_period = active_periods(period_idx, :);
    active_frames = active_period(1):active_period(2);
    plot(time(active_frames), trace(active_frames), 'r');
end
xlim([0 time_window]);
ylim(y_range);
t2 = plot(time(1)*[1 1], y_range, 'k'); % Time indicator
d = plot(time(1), trace(1), 'or',...
            'MarkerFaceColor', 'r',...
            'MarkerSize', 12); % Dot
xlabel('Time [s]');
ylabel('Signal [a.u.]');
hold off;

% Interaction loop:
%   Display the user-specified active period
%------------------------------------------------------------
prompt = 'IC viewer >> ';
resp = lower(strtrim(input(prompt, 's')));
val = str2double(resp);

% State of interaction loop
state.last_val = [];
state.zoomed = false;
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

            otherwise
                fprintf('  Sorry, could not parse "%s"\n', resp);
        end
    end

    resp = lower(strtrim(input(prompt, 's')));
    val = str2double(resp);
end

    % Display subroutine. Note that frames are mean subtracted!
    %------------------------------------------------------------
    function display_active_period(selected_indices)
        for selected_idx = selected_indices
            frames = active_periods(selected_idx,1):...
                     active_periods(selected_idx,2);
            for k = frames
                A = movie(:,:,k);
               % A = A - mean(A(:));
                set(h, 'CData', A);

                % Update time indicators and dot
                set(t1, 'XData', time(k)*[1 1]);
                set(t2, 'XData', time(k)*[1 1]);
                set(d, 'XData', time(k), 'YData', trace(k));

                % Update running trace
                set(a, 'XLim', time(k) + time_window/2*[-1 1]);
                drawnow;
            end
        end
    end % display_active_period

end % main function

function active_frames = parse_active_frames(binary_trace, half_width)
% Segment the active portions of a binary trace into intervals
% 2015 01 31 Tony Hyun Kim

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
end
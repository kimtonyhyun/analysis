function resp = view_ic_over_movie_interactively(ic_filter, time, trace, movie, padding)

% Some parameters
ic_filter_threshold = 0.3;
mad_scale = 8;
active_frame_padding = padding; % Use 100 for 20 Hz movie
time_window = 10; % Width of running window

% Generate the outline of the IC filter
%------------------------------------------------------------
B = threshold_ic_filter(ic_filter, ic_filter_threshold);
B = edge(B, 'canny');
B = ~logical(B);

subplot(3,3,[4 5 7 8]);
% set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]); % Maximize figure
h = imagesc(ic_filter);
colormap gray;
axis image;
xlabel('x [px]');
ylabel('y [px]');

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

% Interaction loop:
%   Display the user-specified active period
%------------------------------------------------------------
prompt = 'IC viewer >> ';
resp = strtrim(input(prompt, 's'));
val = str2double(resp);
while (~isnan(val) || strcmp(resp, 'a')) % Is a number
    if (strcmp(resp, 'a')) % Display all
        display_active_period(1:num_active_periods)
    else
        if ((1 <= val) && (val <= num_active_periods))
            display_active_period(val);
        else
            fprintf('  Sorry, %d is not a valid period index for this IC \n', val);
        end
    end
    resp = strtrim(input(prompt, 's'));
    val = str2double(resp);
end

    % Display subroutine
    %------------------------------------------------------------
    function display_active_period(selected_indices)
        for selected_idx = selected_indices
            frames = active_periods(selected_idx,1):...
                     active_periods(selected_idx,2);
            for k = frames
%                 movie.setDirectory(k);
%                 A = movie.read();
                A = movie(:,:,k);

                % Draw IC edges as black
                A = A - min(A(:));
                A = B.*A;
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
function show_all_traces(ds, varargin)
% Show the normalized activity of all classified cells

% If no cells are classified, then assume all sources are cells
if (ds.num_classified_cells==0)
    cells_to_show = 1:ds.num_cells;
else
    cells_to_show = find(ds.is_cell);
end
num_cells = length(cells_to_show);

amp = 1; % Amplitude of each trace (relative to spacing)
for k = 1:length(varargin)
    vararg = varargin{k};
    if ischar(vararg)
        switch lower(vararg)
            case {'amp', 'amplitude'}
                amp = varargin{k+1};
            case 'cells'
                cells_to_show = varargin{k+1};
        end
    end
end

colors = 'kbr';
num_colors = length(colors);

figure;
hold on;
for k = 1:num_cells
    cell_idx = cells_to_show(k);
    
    tr = ds.get_trace(cell_idx);
    tr_min = min(tr);
    tr_max = max(tr);
    tr = amp*(tr-tr_min)/(tr_max-tr_min);

    color = colors(mod(k,num_colors)+1);
    plot(tr+k-0.5, color);
end
h_ax = gca;

num_frames = length(tr);
y_range = [0 (num_cells-1)+amp];

% DaySummary contains trial structure. In that case, show trial boundaries
%------------------------------------------------------------
if (ds.num_trials > 1)
    trial_start_frames = ds.trial_indices(:,1);
    for k = 1:ds.num_trials
        plot(trial_start_frames(k)*[1 1], y_range, 'k:');
    end
    
    x_ticks = trial_start_frames(1:5:end);
    xticks(x_ticks);
    
    x_tick_labels = 1+5*(0:length(x_ticks)-1);
    xticklabels(x_tick_labels);
    xlabel(sprintf('Trials (%d total)', ds.num_trials));
else
    xlabel(sprintf('Frames (%d total)', num_frames));
end
hold off;

xlim([1 num_frames]);
yticks(1:num_cells);
ylabel(sprintf('Cell index (%d total)', num_cells));
set(gca, 'TickLength', [0 0]);

% Navigation controls
% We use "UserData" to encode GUI state
%------------------------------------------------------------
zoom_all();

bw = 0.06; % Button width
bh = 1/6;   % Button height

uicontrol('Style', 'pushbutton',...
    'String', 'Down',...
    'Units', 'normalized',...
    'Position', [1-bw 0 bw bh],...
    'Callback', @scroll_down);

uicontrol('Style', 'pushbutton',...
    'String', '5',...
    'Units', 'normalized',...
    'Position', [1-bw bh bw bh],...
    'Callback', @zoom5);

uicontrol('Style', 'pushbutton',...
    'String', '10',...
    'Units', 'normalized',...
    'Position', [1-bw 2*bh bw bh],...
    'Callback', @zoom10);

uicontrol('Style', 'pushbutton',...
    'String', '25',...
    'Units', 'normalized',...
    'Position', [1-bw 3*bh bw bh],...
    'Callback', @zoom25);

uicontrol('Style', 'pushbutton',...
    'String', 'All',...
    'Units', 'normalized',...
    'Position', [1-bw 4*bh bw bh],...
    'Callback', @zoom_all);

uicontrol('Style', 'pushbutton',...
    'String', 'Up',...
    'Units', 'normalized',...
    'Position', [1-bw 5*bh bw bh],...
    'Callback', @scroll_up);

    function zoom5(~,~)
        current_range = h_ax.UserData;
        h_ax.UserData = current_range(1) + [0 4];
        zoom_to_cells();
    end

    function zoom10(~, ~)
        current_range = h_ax.UserData;
        h_ax.UserData = current_range(1) + [0 9];
        zoom_to_cells();
    end

    function zoom25(~, ~)
        current_range = h_ax.UserData;
        h_ax.UserData = current_range(1) + [0 24];
        zoom_to_cells();
    end

    function zoom_all(~, ~)
        h_ax.UserData = [1 num_cells];
        zoom_to_cells();
    end

    function scroll_down(~, ~)
        current_range = h_ax.UserData;
        if current_range(1) > 1
            h_ax.UserData = current_range - (diff(current_range) + 1);
            zoom_to_cells();
        end
    end

    function scroll_up(~, ~)
        current_range = h_ax.UserData;
        if current_range(2) < num_cells
            h_ax.UserData = current_range + (diff(current_range) + 1);
            zoom_to_cells();
        end
    end

    function zoom_to_cells()
        ylim(h_ax, h_ax.UserData + [-0.5 0.5]);
    end

end % show_all_traces
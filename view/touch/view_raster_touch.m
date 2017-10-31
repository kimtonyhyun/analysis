function view_raster_touch(ds, cell_idx, varargin)
% Tool for browsing single cell rasters of a single day (i.e. DaySummary),
% but without keyboard interaction!

h_fig = [];
subraster_type = '';

for i = 1:length(varargin)
    vararg = varargin{i};
    if ischar(vararg)
        switch lower(vararg)
            case 'fig'
                h_fig = varargin{i+1};
            case 'type'
                subraster_type = varargin{i+1};
        end
    end
end

if isempty(h_fig)
    h_fig = figure;
else
    figure(h_fig);
    clf;
end

% Main drawing routine
%------------------------------------------------------------
% Image of cell
subplot(3,2,1);
imagesc(ds.cells(cell_idx).im);
axis image;
title(sprintf('Cell %d (%s)', cell_idx, ds.cells(cell_idx).label));
colormap jet; freezeColors;

% Raster of all trials, with correctness
h_full_raster = subplot(3,2,[3 5]);
ds.plot_cell_raster(cell_idx, 'draw_correct');
raster_scale = get(gca, 'CLim'); % Scale that applies to all trials
title('All trials');
colormap jet; freezeColors;
if ds.is_switchdata_loaded
    x_ends = get(gca, 'XLim');
    hold on;
    plot(x_ends, (ds.switchdata.pre_switch_trials(end)+0.5)*[1 1], 'w--', 'LineWidth', 1);
    plot(x_ends, (ds.switchdata.post_switch_trials(1)-0.5)*[1 1], 'w--', 'LineWidth', 1);
    hold off;
end
set(h_full_raster, 'ButtonDownFcn', @select_trial);

% Default subraster option
if isempty(subraster_type)
    if ds.is_switchdata_loaded
        subraster_type = 'path';
    else
        subraster_type = 'standard';
    end
end
draw_subraster(subraster_type);

% Navigation controls
%------------------------------------------------------------
back_btn = uicontrol('Style', 'pushbutton',...
    'String', '<<',...
    'Units', 'normalized',...
    'Position', [0 0.1 0.05 0.8],...
    'Callback', @back_to_browse);

% NOTE: The up-down navigation assumes that the DaySummary contains no
% non-cells
up_btn = uicontrol('Style', 'pushbutton',...
    'String', '^^',...
    'Units', 'normalized',...
    'Position', [0 0.9 0.05 0.1],...
    'Callback', {@jump_to_cell, cell_idx-1});
if (cell_idx == 1)
    down_btn.Enable = 'off';
end
down_btn = uicontrol('Style', 'pushbutton',...
    'String', 'vv',...
    'Units', 'normalized',...
    'Position', [0 0 0.05 0.1],...
    'Callback', {@jump_to_cell, cell_idx+1});
if (cell_idx == ds.num_cells)
    down_btn.Enable = 'off';
end

% Visualize alternate rasters
std_raster_btn = uicontrol('Style', 'pushbutton',...
    'String', 'S',...
    'TooltipString', 'Standard subrasters',...
    'Units', 'normalized',...
    'Position', [0.95 0.9 0.05 0.1],...
    'Callback', {@redraw_callback, 'standard'});
path_raster_btn = uicontrol('Style', 'pushbutton',...
    'String', 'P',...
    'TooltipString', 'Path-specific subrasters',...
    'Units', 'normalized',...
    'Position', [0.95 0.8 0.05 0.1],...
    'Callback', {@redraw_callback, 'path'});
constant_path_btn = uicontrol('Style', 'pushbutton',...
    'String', 'CoP',...
    'TooltipString', 'Constant path analysis',...
    'Units', 'normalized',...
    'Position', [0.95 0.7 0.05 0.1],...
    'Callback', {@redraw_callback, 'constant_path'});
changing_path_btn = uicontrol('Style', 'pushbutton',...
    'String', 'ChP',...
    'TooltipString', 'Changing path analysis',...
    'Units', 'normalized',...
    'Position', [0.95 0.6 0.05 0.1],...
    'Callback', {@redraw_callback, 'changing_path'});
if ~ds.is_switchdata_loaded
    path_raster_btn.Enable = 'off';
    constant_path_btn.Enable = 'off';
    changing_path_btn.Enable = 'off';
end

    function back_to_browse(~, ~)
        browse_rasters_touch(ds, 'fig', h_fig, 'cell_idx', cell_idx);
    end

    function jump_to_cell(~, ~, new_idx)
        % Maintain subraster type when browsing cells
        view_raster_touch(ds, new_idx, 'fig', h_fig, 'type', subraster_type);
    end

    function select_trial(~, e)
        trial_idx = round(e.IntersectionPoint(2));
        trial_idx = max(1, trial_idx);
        trial_idx = min(trial_idx, ds.num_trials);
        if ds.is_behavior_loaded
            view_detailed_trial_touch(ds, cell_idx, trial_idx, 'fig', h_fig);
        else
            fprintf('Error: Behavior video has not been loaded into this DaySummary!\n');
        end
    end

    function redraw_callback(~, ~, type)
        draw_subraster(type);
    end

    function draw_subraster(type)
        switch type
            case 'standard'
                draw_standard_subrasters(ds, cell_idx, raster_scale);
                subraster_type = type;
            case 'path'
                draw_path_subrasters(ds, cell_idx, raster_scale);
                subraster_type = type;
            case 'constant_path'
                draw_constant_path_analysis(ds, cell_idx);
                subraster_type = type;
            otherwise
                msgbox(sprintf('Type "%s" not implemented!', type));
        end
    end

end % draw_raster_page

function view_raster_touch(ds, cell_idx, varargin)
% Tool for browsing single cell rasters of a single day (i.e. DaySummary),
% but without keyboard interaction!

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
set(h_full_raster, 'ButtonDownFcn', @select_trial);

% Divide rasters by correctness
% subplot(3,4,3);
% ds.plot_cell_raster(cell_idx, 'correct');
% set(gca, 'CLim', raster_scale);
% title('Correct');
% subplot(3,4,4);
% ds.plot_cell_raster(cell_idx, 'incorrect');
% set(gca, 'CLim', raster_scale);
% title('Incorrect');

% Divide rasters by start location
% subplot(3,4,7);
% ds.plot_cell_raster(cell_idx, 'start', 'west', 'draw_correct');
% set(gca, 'CLim', raster_scale);
% title('West start'); 
% subplot(3,4,8);
% ds.plot_cell_raster(cell_idx, 'start', 'east', 'draw_correct');
% set(gca, 'CLim', raster_scale);
% title('East start');
% 
% % Divide rasters by end location
% subplot(3,4,11);
% ds.plot_cell_raster(cell_idx, 'end', 'south', 'draw_correct');
% set(gca, 'CLim', raster_scale);
% title('South end');
% subplot(3,4,12);
% ds.plot_cell_raster(cell_idx, 'end', 'north', 'draw_correct');
% set(gca, 'CLim', raster_scale);
% title('North end');

% Show constant path rasters
subplot(1,4,3);
ds.plot_cell_raster(cell_idx, 'start', 'east', 'end', 'north', 'draw_correct');
set(gca, 'CLim', raster_scale);
title('Constant path: East-North');

% Show changing path rasters
subplot(2,4,4);
ds.plot_cell_raster(cell_idx, 'start', 'west', 'end', 'south', 'draw_correct');
set(gca, 'CLim', raster_scale);
title('Changing path - PRE: West-South');

subplot(2,4,8);
ds.plot_cell_raster(cell_idx, 'start', 'west', 'end', 'north', 'draw_correct');
set(gca, 'CLim', raster_scale);
title('Changing path - POST: West-North');

% Navigation controls
%------------------------------------------------------------
back_btn = uicontrol('Style', 'pushbutton',...
    'String', '<<',...
    'Units', 'normalized',...
    'Position', [0 0 0.05 1],...
    'Callback', @back_to_browse);

    function back_to_browse(~, ~)
        browse_rasters_touch(ds, 'fig', h_fig, 'cell_idx', cell_idx);
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

end % draw_raster_page

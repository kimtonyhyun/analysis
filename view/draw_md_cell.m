function draw_md_cell(md, common_cell_idx, raster_fn)

% Original behavior: Use the built-in raster visualization function of DS
if ~exist('raster_fn', 'var')
    raster_fn = @(ds,k) plot_cell_raster(ds, k, 'draw_correct');
end

% We allow 'raster_fn' to be a cell array of multiple raster visualization
% functions. The following allows for a common call interface, whether 1 or
% more raster_fns are provided.
%
% Each 'raster_fn' should be a function handle with the arguments @(ds,k)
% where ds is the DaySummary and k is the cell index.
if length(raster_fn) == 1
    raster_fn = {raster_fn};
end
num_rasters = length(raster_fn);

if num_rasters > 1
    sp_num_rows = 1 + 2*num_rasters;
else % num_rasters == 1
    sp_num_rows = 4;
end

sp = @(m,n,p) subtightplot(m, n, p, 0.05, 0.05, 0.05); % Gap, Margin-X, Margin-Y

traces = cell(1, md.num_days);
for k = 1:md.num_days
    day = md.valid_days(k);
    day_name = md.valid_day_names{k};
    cell_idx_k = md.get_cell_idx(common_cell_idx, day);
    ds = md.day(day); % Just a shorthand

    % Draw image of cell on each day
    sp(sp_num_rows, md.num_days, k);
    cell_image_k = ds.cells(cell_idx_k).im;
    cell_com_k = ds.cells(cell_idx_k).com;
    zoom_half_width = min(size(cell_image_k))/10;

    imagesc(cell_image_k);
    axis image;
    xlim(cell_com_k(1)+zoom_half_width*[-1 1]);
    ylim(cell_com_k(2)+zoom_half_width*[-1 1]);

    title_str = sprintf('%s\nCell %d', day_name, cell_idx_k);
    if (day == md.sort_day)
        title(title_str, 'FontWeight', 'bold', 'FontSize', 12);
    else
        title(title_str, 'FontWeight', 'normal');
    end

    % Draw raster
    for r = 1:num_rasters
        sp(sp_num_rows, md.num_days,...
            [(1+2*(r-1))*md.num_days+k (2+2*(r-1))*md.num_days+k]);
        raster_fn{r}(ds, cell_idx_k);
    end

    % Get trace
    traces{k} = ds.get_trace(cell_idx_k);
end

% If we're plotting only 1 raster per cell, then plot cross-day traces as a
% single plot in the bottom row. If there are more than 1 raster per cell,
% we skip this in order to conserve vertical space
if num_rasters == 1
    trace_offsets = zeros(1+md.num_days, 1);
    for k = 1:md.num_days
        trace_offsets(k+1) = length(traces{k});
    end
    trace_offsets = cumsum(trace_offsets);
    
    colors = 'rbk';
    sp(sp_num_rows, md.num_days,...
        ((sp_num_rows-1)*md.num_days+1):(sp_num_rows*md.num_days));
    for k = 1:md.num_days
        color = colors(1+mod(k,length(colors)));
        plot(trace_offsets(k):(trace_offsets(k+1)-1),...
             traces{k}, color);
        hold on;
    end
    hold off;
    grid on;
    xlim([0 trace_offsets(end)-1]);
    set(gca, 'XTick', []);
    full_trace = cell2mat(traces);
    M = max(full_trace);
    m = min(full_trace);
    ylim([m M] + 0.1*(M-m)*[-1 1]);
    ylabel('Trace [a.u.]');
end
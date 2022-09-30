function draw_md_cell(md, common_cell_idx)

sp = @(m,n,p) subtightplot(m, n, p, 0.05, 0.05, 0.05); % Gap, Margin-X, Margin-Y

traces = cell(1, md.num_days);
for k = 1:md.num_days
    day = md.valid_days(k);
    cell_idx_k = md.get_cell_idx(common_cell_idx, day);
    ds = md.day(day); % Just a shorthand

    % Draw image of cell on each day
    sp(4, md.num_days, k);
    cell_image_k = ds.cells(cell_idx_k).im;
    cell_com_k = ds.cells(cell_idx_k).com;
    zoom_half_width = min(size(cell_image_k))/10;

    imagesc(cell_image_k);
    axis image;
    xlim(cell_com_k(1)+zoom_half_width*[-1 1]);
    ylim(cell_com_k(2)+zoom_half_width*[-1 1]);

    title_str = sprintf('Day %d -- Cell %d', day, cell_idx_k);
    if (day == md.sort_day)
        title(title_str, 'FontWeight', 'bold', 'FontSize', 12);
    else
        title(title_str, 'FontWeight', 'normal');
    end

    % Draw raster
    sp(4, md.num_days, [md.num_days+k 2*md.num_days+k]);
    % Programmatically specify the raster visualization function
    ctxstr.vis.plot_raster_from_ds(ds, cell_idx_k, 2);
%             ds.plot_cell_raster(cell_idx_k, 'draw_correct');

    % Get trace
    traces{k} = ds.get_trace(cell_idx_k);
end

% Draw traces on a single plot
trace_offsets = zeros(1+md.num_days, 1);
for k = 1:md.num_days
    trace_offsets(k+1) = length(traces{k});
end
trace_offsets = cumsum(trace_offsets);

colors = 'rbk';
sp(4, md.num_days, (3*md.num_days+1):(4*md.num_days));
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

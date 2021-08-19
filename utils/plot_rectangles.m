function plot_rectangles(x_ranges, y_lims)

for k = 1:size(x_ranges,1)
    x_range = x_ranges(k,:);
    w = x_range(2) - x_range(1);
    h = y_lims(2) - y_lims(1);
    rectangle('Position', [x_range(1) y_lims(1) w h],...
        'EdgeColor', 'none', 'FaceColor', [0 1 1 0.15],...
        'HitTest', 'off');
end

function plot_depth_boundaries(ds_list)
% Plot boundaries of DaySummaries in 'ds_list' (cell of DaySummaries)

num_depths = length(ds_list);
color_offset = 0; % Avoid bright colors of the colormap
colors = parula(num_depths+color_offset);
colors = flipud(colors(1:num_depths,:));

plot_boundaries_with_transform(ds_list{1}, colors(1,:));
hold on;
for k = 2:num_depths
    plot_boundaries_with_transform(ds_list{k}, colors(k,:));
end
hold off;
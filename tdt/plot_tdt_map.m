function plot_tdt_map(ds)
% Wrapper for 'plot_cell_map' for displaying tdTomato identities

[pos, neg] = collect_tdt_labels(ds);
unlabeled = setdiff(1:ds.num_cells, [pos neg]);

color_grouping = {pos, [0.85 0.325 0.098]; % Orange
                  neg, [0 0.447 0.741]; % Blue
                  unlabeled, 'w'};

ds.plot_cell_map(color_grouping);

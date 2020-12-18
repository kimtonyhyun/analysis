function plot_cell_map_redblue(obj, cell_inds, cell_values, varargin)
% Plot cell map, but with each cell's color derived from an analog value.
% Conversion from value to color uses the redblue colormap, where the 
% midpoint of the redblue colormap (white) maps to a value of zero.
%
% Example use cases:
%   - Color the cell map based on correlation scores
%   - Color the cell map based on the number of detected events

% We use 100 shades of colors on each side of 0
colors = redblue(201);

max_value = max(abs(cell_values));
cell_values = cell_values / max_value; % Rescale to [-1, 1]

num_cells = length(cell_inds);
color_grouping = cell(num_cells, 2);
for k = 1:num_cells
    cell_ind = cell_inds(k);
    cell_value = cell_values(k);
    
    color_ind = round(cell_value*100)+101;
    
    color_grouping{k,1} = cell_ind;
    color_grouping{k,2} = colors(color_ind,:);
end

obj.plot_cell_map(color_grouping);

end
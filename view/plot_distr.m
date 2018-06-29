function [X, Y] = plot_distr(ys, color, x_offset)

num_groups = length(ys);

% Compute the 25th, 50th (median), 75th percentiles of each group
quartiles = cellfun(@(y) prctile(y,[25 50 75]), ys, 'UniformOutput', false);
quartiles = cell2mat(quartiles);

% Plot the 25-75 range
x = 1:num_groups;
X = kron(x, [1 1 NaN]);
Y = [quartiles(:,[1 3]) NaN(num_groups,1)]';
Y = Y(:);
plot(X + x_offset, Y, 'Color', color);
hold on;

% Plot the median
plot(x + x_offset, quartiles(:,2), '.', 'Color', color);

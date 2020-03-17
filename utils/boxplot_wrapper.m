function [d, g] = boxplot_wrapper(grouped_data, varargin)
% A wrapper for 'boxplot' that allows a cell-based interface, where:
%   - grouped_data: [num_groups x 2] cell, where
%       - grouped_data{k,1}: Name of group. All numeric or all string.
%       - grouped_data{k,2}: Data belonging to group. Needs to be column
%       vector

num_groups = size(grouped_data, 1);
num_data_per_group = zeros(num_groups, 1);

for k = 1:num_groups
    num_data_per_group(k) = length(grouped_data{k,2});
end
num_all_data = sum(num_data_per_group);

% Interface to boxplot
d = cell2mat(grouped_data(:,2));
g = cell(num_all_data, 1);
idx = 1;
for k = 1:num_groups
    g(idx:idx+num_data_per_group(k)-1) = grouped_data(k,1);
    idx = idx + num_data_per_group(k);
end

boxplot(d, g, varargin{:});
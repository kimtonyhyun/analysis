function [d, g] = boxplot_wrapper(group_names, group_data, varargin)
% A wrapper for 'boxplot' that allows a cell-based interface, where:
%   - group_names: Names of groups. All numeric or all string.
%   - group_data: Data belonging to each group.
%
% Note that 'boxplot' will only sort the groups in increasing group_name
% if the names are numeric!

group_names = group_names(:); % Force column vector
group_data = group_data(:);

num_groups = size(group_names, 1);
num_data_per_group = zeros(num_groups, 1);

for k = 1:num_groups
    gdk = group_data{k};
    num_data_per_group(k) = length(gdk);
    
    % Force column vector. Note that both the containing cell and its 
    % contents are column vectors. Needed for 'cell2mat', below.
    group_data{k} = gdk(:);
end
num_all_data = sum(num_data_per_group);

% Interface to boxplot
d = cell2mat(group_data);
g = cell(num_all_data, 1);
idx = 1;
for k = 1:num_groups
    if iscell(group_names)
        g(idx:idx+num_data_per_group(k)-1) = group_names(k);
    else
        g(idx:idx+num_data_per_group(k)-1) = {group_names(k)};
    end
    idx = idx + num_data_per_group(k);
end

boxplot(d, g, varargin{:});
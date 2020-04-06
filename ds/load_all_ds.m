function [ds_list, num_ds] = load_all_ds(dataset_stem, subpath_to_rec, get_id_fn)
% Convenience function for loading multiple DaySummary's.
%
% Inputs:
%   - 'dataset_stem': Each subdirectory containing a DaySummary must start
%                     with 'dataset_stem'.
%   - 'subpath_to_rec': Subpath to the rec file.
%   - 'get_id_fn': Function (string --> string) to extract the "id" of
%                each DaySummary

path_to_ds = fileparts(dataset_stem);

datasets = dir(sprintf('%s*', dataset_stem));
num_ds = length(datasets);

ds_list = cell(num_ds, 2); % [ID(str) DaySummary]
for k = 1:num_ds
    dataset_name = datasets(k).name;
    path_to_rec = fullfile(path_to_ds, dataset_name, subpath_to_rec);
    
    ds_list{k,1} = get_id_fn(dataset_name);
    ds_list{k,2} = DaySummary([], path_to_rec);
end
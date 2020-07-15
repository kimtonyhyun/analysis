function [raster, info] = compute_averaged_raster(ds)
% Generate raster data where each row is the trial-averaged fluorescence
% trace of each cell.
%
% Can view the results as:
%   imagesc(info.aligned_time, 1:info.num_cells, raster);
%

cell_inds = find(ds.is_cell);
num_cells = length(cell_inds);

% FIXME: Allow for selection of trials
trial_inds = find([ds.trials.correct]);

% FIXME: Allow for alignment to other trial timepoints
align_to = 3;
alignment_frames = ds.trial_indices(trial_inds, align_to);

[~, r_info] = ds.get_aligned_trace(1, trial_inds, alignment_frames);
num_frames_per_trial = length(r_info.aligned_time);

raster = zeros(num_cells, num_frames_per_trial);
max_frame_inds = zeros(num_cells, 1);

for k = 1:num_cells
    cell_ind = cell_inds(k);
    
    % First, get the cell raster: [num_trials x num_frames_per_trial]
    r = ds.get_aligned_trace(cell_ind, trial_inds, alignment_frames);
    avg_tr = mean(r,1);
    
    % Next, normalize the average trace
    M = max(avg_tr);
    m = min(avg_tr);
    avg_tr = (avg_tr-m)/(M-m);
    raster(k,:) = avg_tr;
    
    % Finally, compute the frame idx where the avg trace is maximal (for
    % sorting of the raster)
    [~, max_frame_ind] = max(avg_tr);
    max_frame_inds(k) = max_frame_ind;
end

% Sort rows by time to fluorescence peak
[~, sort_ind] = sort(max_frame_inds, 'ascend');
raster = raster(sort_ind, :);

% Package for output
info.orig_cell_inds = sort_ind;
info.num_cells = num_cells;
info.aligned_time = r_info.aligned_time;
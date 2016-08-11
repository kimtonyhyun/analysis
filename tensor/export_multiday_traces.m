function [X, neuron_map, trial_map] = export_multiday_traces(md, trial_type)
% Exports cross-day aligned cell traces (via MultiDay) into a lightweight
% format (description below) for tensor analysis
%
% Output format:
%   X: Column cell vector [M x 1] where M is the total number of trials
%      (across all days) of the requested 'trial_type' (e.g. 'en')
%       
%      X{j} is a matrix of cell traces, in the format [num_matched_cells x
%      num_frames_in_trial_j].
%
%   neuron_map: Matrix [num_matched_cells x num_days] that maps each 
%      matched neuron to its per-day neuron index.
%
%   trial_map: Matrix [M x 2] where M is the total number of trials of 
%      requested 'trial_type'. For the i-th trial, K(i,1) is the day index 
%      of that trial, and K(i,2) is the trial index in the original day.
%   
max_num_trials = 0;
for day_idx = md.valid_days
    max_num_trials = max_num_trials + md.day(day_idx).num_trials;
end

neuron_map = md.matched_indices;

% Preallocate outputs
X = cell(max_num_trials,1);
trial_map = zeros(max_num_trials, 2); % Format: [Day-index, Trial-index]

idx = 0;
for k = 1:md.num_days
    day_idx = md.valid_days(k);
    day = md.day(day_idx);
    day_cell_indices = md.matched_indices(:,k);
    
    % Apply trial type filter. TODO: Factor out the same bit of code in
    % 'session_tensor' to make use of the same helper function
    filtered_trials = filter_trials_tensor(day, trial_type);
    
    for trial_idx = filtered_trials
        traces = day.trials(trial_idx).traces; % [neuron x time]
        
        idx = idx + 1;
        X{idx} = traces(day_cell_indices,:); % Reordered to use common index
        trial_map(idx,:) = [day_idx, trial_idx];
    end
end

% Truncate outputs
X = X(1:idx);
trial_map = trial_map(1:idx,:);
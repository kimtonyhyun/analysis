function X = export_multiday_traces(md)

max_num_trials = 0;
for day_idx = md.valid_days
    max_num_trials = max_num_trials + md.day(day_idx).num_trials;
end

X = cell(max_num_trials,1); % Preallocate output

idx = 0;
for k = 1:md.num_days
    day_idx = md.valid_days(k);
    day = md.day(day_idx);
    cell_indices = md.matched_indices(:,k);
    
    for trial_idx = 1:day.num_trials
        traces = day.trials(trial_idx).traces; % [neuron x time]
        
        idx = idx + 1;
        X{idx} = traces(cell_indices,:); % Reordered to use common index
    end
end
function events_per_trial = compute_eventdata_per_trial(eventdata, trial_indices)

num_trials = size(trial_indices, 1);
events_per_trial = cell(num_trials, 1);

frames2trial = zeros(trial_indices(end,end), 1);
for k = 1:num_trials
    frames2trial(trial_indices(k,1):trial_indices(k,end)) = k;
end

num_events = size(eventdata,1);
events2trial = frames2trial(eventdata(:,2));
new_trial_inds = [1; find(diff(events2trial))+1];
n = length(new_trial_inds);

for k = 1:n
    start_idx = new_trial_inds(k);
    if k ~= n
        end_idx = new_trial_inds(k+1)-1;
    else
        end_idx = num_events;
    end
    
    tidx = events2trial(start_idx);
    events_per_trial{tidx} = eventdata(start_idx:end_idx, :);
end

end % compute_eventdata_per_trial
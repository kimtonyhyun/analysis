function events = find_events_in_trials(trace, trial_indices, threshold, baseline)

num_trials = size(trial_indices, 1);

events = cell(num_trials, 1);
for k = 1:num_trials
    % Find events within the trial
    trial_frames = trial_indices(k,[1 end]);
    trial_frames = trial_frames(1):trial_frames(2);
    es = find_events(trace(trial_frames), threshold, baseline);
    
    % Add frame offset for the k-th trial
    if ~isempty(es)
        es(:,1) = es(:,1) + trial_frames(1) - 1;
        es(:,2) = es(:,2) + trial_frames(1) - 1;
    end
    events{k} = es;
end

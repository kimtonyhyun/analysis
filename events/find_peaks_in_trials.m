function peaks = find_peaks_in_trials(trace, trial_indices, threshold, baseline)
% Detect peaks in the fluorescence trace. Wrapper for 'find_peaks' that is
% aware of trial boundaries.

num_trials = size(trial_indices, 1);

peaks = cell(num_trials, 1);
for k = 1:num_trials
    % Find all fluorescence peaks within each trial
    trial_frames = trial_indices(k,[1 end]);
    trial_frames = trial_frames(1):trial_frames(2);
    peaks_k = find_peaks(trace(trial_frames), threshold, baseline);
    
    % Add frame offset for the k-th trial
    frame_offset = trial_frames(1) - 1;
    if ~isempty(peaks_k)
        peaks_k(:,1) = peaks_k(:,1) + frame_offset; % Pre-trough frame
        peaks_k(:,2) = peaks_k(:,2) + frame_offset; % Peak frame
    end
    peaks{k} = peaks_k;
end
peaks = cell2mat(peaks);

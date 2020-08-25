function events = find_events_in_trials(trace, trial_indices, threshold, baseline, amp_threshold)
% Detect events in the fluorescence trace. Expects a smoothed version of
% the fluorescence trace.
%
% Returns:
%   events: [num_events x 3] where
%       events(k,1): Frame of the trough preceding the k-th event.
%       events(k,2): Frame of the peak of the k-th event.
%       events(k,3): Peak height relative to preceding trough ("Event amplitude")
%
% Basic usage:
%   trace = filter_trace(ds, cell_idx, fps, cutoff_freq);
%   [baseline, sigma] = estimate_baseline_sigma(trace);
%   events = find_events_in_trials(trace, ds.trial_indices, baseline+5*sigma, baseline, 0.1);
%
% Can also be used interactively:
%   events = detect_events(ds, cell_idx, fps);
%

num_trials = size(trial_indices, 1);

% First, we find all peaks (local maxima) in the fluorescence trace that
% exceed the given threshold.
%------------------------------------------------------------
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
        peaks_k(:,3) = peaks_k(:,3) + frame_offset; % Peak frame
        peaks_k(:,5) = peaks_k(:,5) + frame_offset; % Post-trough frame
    end
    peaks{k} = peaks_k;
end
peaks = cell2mat(peaks);

if isempty(peaks)
    events = [];
    return;
end

% Output format
% TODO: Consider named fields (e.g. struct)?
events(:,1) = peaks(:,1); % Pre-trough frame
events(:,2) = peaks(:,3); % Peak frame
events(:,3) = peaks(:,4) - peaks(:,2); % Peak height relative to pre-trough

% Filter for amplitude heights
if ~isempty(events)
    max_event_amplitude = max(events(:,3));
    filtered_events = events(:,3) > amp_threshold * max_event_amplitude;
    events = events(filtered_events,:);
end
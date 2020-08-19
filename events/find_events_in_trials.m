function events = find_events_in_trials(trace, trial_indices, threshold, baseline, merge_threshold, amp_threshold)
% Detect events in the fluorescence trace. Expects a smoothed version of
% the fluorescence trace.
%
% Returns:
%   events: [num_events x 3] where
%       events(k,1): Frame of the trough preceding the k-th event.
%       events(k,2): Frame of the peak of the k-th event.
%       events(k,3): Peak height relative to preceding trough ("Event amplitude")
%       events(k,4): Peak height relative to subsequent trough
%
% Basic usage:
%   trace = filter_trace(ds, cell_idx, fps, cutoff_freq);
%   [baseline, sigma] = estimate_baseline_sigma(trace);
%   events = find_events_in_trials(trace, ds.trial_indices, baseline+5*sigma, baseline, 0.1, 0.1);
%
% Can also be used interactively:
%   events = detect_events(ds, cell_idx, fps);
%

if (nargin < 5)
    amp_threshold = 0;
end

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
num_peaks = size(peaks, 1);

if isempty(peaks)
    events = [];
    return;
end

% Next, process the list of all peaks
%------------------------------------------------------------
max_trace_amplitude = max(peaks(:,4)) - baseline;

% Merge nearby and "indistinct" peaks. See below for the exact criteria.
merged_peaks = zeros(num_peaks, 6);
idx = 1;
peak_to_be_written = peaks(1,:);
for k = 2:num_peaks
    p = peaks(k,:);
    % First, the peaks need to "connect", i.e. the POST trough frame of the
    % preceding peak is equal to the PRE trough frame of the next peak.
    merge_peak = false;
    if peak_to_be_written(5) == p(1)
        % Secondly, successive peak heights mut be increasing
        if peak_to_be_written(4) < p(4)
            % Thirdly, successive pre-trough values must be increasing
            if peak_to_be_written(2) < p(2)
                % Finally, the POST trough amplitude must be smaller than the
                % merge threshold.
                post_trough_amp = peak_to_be_written(4) - peak_to_be_written(6);
                if post_trough_amp < merge_threshold * max_trace_amplitude
                    merge_peak = true;
                end
            end
        end
    end
    
    if merge_peak
        % Keep PRE trough values, but replace peak and POST trough values
        peak_to_be_written(3:6) = p(3:6);
    else
        merged_peaks(idx,:) = peak_to_be_written;
        idx = idx + 1;
        
        peak_to_be_written = p;
    end
end
merged_peaks(idx,:) = peak_to_be_written;
merged_peaks = merged_peaks(1:idx,:);

% Output format
% TODO: Consider named fields (e.g. struct)?
events(:,1) = merged_peaks(:,1); % Pre-trough frame
events(:,2) = merged_peaks(:,3); % Peak frame
events(:,3) = merged_peaks(:,4) - merged_peaks(:,2); % Peak height relative to pre-trough
events(:,4) = merged_peaks(:,4) - merged_peaks(:,6); % Peak height relative to post-trough

% Filter for amplitude heights
if ~isempty(events)
    max_event_amplitude = max(events(:,3));
    filtered_events = events(:,3) > amp_threshold * max_event_amplitude;
    events = events(filtered_events,:);
end
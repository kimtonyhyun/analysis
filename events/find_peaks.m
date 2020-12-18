function peaks = find_peaks(trace, threshold, baseline)
% For each segment in 'trace' above the specified 'threshold', every local
% maximum within those segments are identified as peaks.
%
% Returns:
%   peaks: [num_peaks x 3] where
%     peaks(k,1): Frame of the trough preceding the k-th peak. ("Pre" trough)
%     peaks(k,2): Frame of the k-th peak.
%     peaks(k,3): Amplitude of the k-th peak (peak - trough)
%

above_thresh_frames = find(trace >= threshold);

if ~any(above_thresh_frames)
    peaks = [];
    return;
end

segments = frame_list_to_segments(above_thresh_frames);
num_segments = size(segments,1);

peak_frames = cell(1,num_segments);
for k = 1:num_segments
    seg = segments(k,1):segments(k,2);
    tr_seg = trace(seg);
    
    peak_frames{k} = find_all_localmax(seg, tr_seg);
end
peak_frames = cell2mat(peak_frames);

num_peaks = length(peak_frames);
peaks = zeros(num_peaks, 3);
for k = 1:num_peaks
    peak_frame = peak_frames(k);
    
    % PRE trough
    %------------------------------------------------------------
    % Attempt to find the trough preceding the event peak. Note that there
    % are cases where the trough cannot be found:
    %   1) The peak is at the beginning of the trace
    %   2) The trough is at the beginning of the trace -- which means that
    %   it's possible that we didn't roll all the way down to the trough
    pre_trough_not_found = false;
    if (peak_frame == 1)
        pre_trough_not_found = true;
    else
        pre_trough_frame = seek_localmin(trace, peak_frame - 1);
        if (pre_trough_frame == 1)
            pre_trough_not_found = true;
        end
    end
    
    % In the case that trough was not found, use the baseline
    if pre_trough_not_found
        pre_trough_frame = 1;
        pre_trough_val = baseline;
    else
        pre_trough_val = trace(pre_trough_frame); 
    end
    
    % Format data for output
    %------------------------------------------------------------
    peaks(k,1) = pre_trough_frame;
    peaks(k,2) = peak_frame;
    peaks(k,3) = trace(peak_frame) - pre_trough_val;
end

end % find_events

function localmax = find_all_localmax(x_seg, y_seg)
    % A point is a local maximum if it is larger than the value to the left
    % and to the value to the right.
    
    % Whether the first and last points of the segment can be
    % identified as a local maximum
    allow_ends = true;
    x1 = [allow_ends y_seg(2:end)>y_seg(1:end-1)];
    x2 = [y_seg(1:end-1)>y_seg(2:end) allow_ends];

    localmax = x_seg(x1 & x2);
end % find_all_localmax

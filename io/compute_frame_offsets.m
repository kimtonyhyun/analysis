function [pre_offset, post_offset] = compute_frame_offsets(frame_indices, align_index)
    % Compute the maximum length of pre- and post-alignment frames that is
    % common to all trials
    frame_indices = double(frame_indices);
    
    alignment_frames = frame_indices(:,align_index);
    frame_indices = frame_indices - repmat(alignment_frames, [1 4]); % Align each trial to closing of gate
    
    pre_offset = max(frame_indices(:,1));
    post_offset = min(frame_indices(:,4));
end
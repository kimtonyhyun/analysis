function [ frame_indices ] = make_frame_indices( num_frames )
% Creates frame indices so that get_mouse_XY_pos can go
%   through the behavior video in chunks (not enough memory to
%   load the whole video)
%
%   splits up frames into chunks of 500
%   returns a matrix like [1    500
%                          501  1000
%                          1001 1500 ...] etc. until the last frame
%
% 2015-02-26 Fori Wang

    frame_indices = [];
    for chunk_idx = 1:500:num_frames
        if chunk_idx+499 < num_frames
            frame_indices = [frame_indices; chunk_idx chunk_idx+499];
        else % for the last chunk        
            frame_indices = [frame_indices; chunk_idx num_frames];
        end
    end

end


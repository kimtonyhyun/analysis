function [ frame_indices ] = make_frame_indices( num_frames , chunk_size)
% Creates frame indices so that get_mouse_XY_pos can go
%   through the behavior video in chunks (not enough memory to
%   load the whole video)
%
%   splits up frames into chunks of size chunk_size
%     For example, [frame_indices] = make_frame_indices(50000, 1000)
%       returns:           [1    1000
%                          1001  2000
%                          2001  3000 ...] etc. until the last frame
%
% 2015-02-26 Fori Wang

    frame_indices = [];
    for chunk_idx = 1:chunk_size:num_frames
        if chunk_idx+chunk_size-1 < num_frames
            frame_indices = [frame_indices; chunk_idx chunk_idx+chunk_size-1];
        else % for the last chunk        
            frame_indices = [frame_indices; chunk_idx num_frames];
        end
    end

end


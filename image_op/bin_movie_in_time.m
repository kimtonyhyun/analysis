function M_b = bin_movie_in_time(M, bin_factor)
% Temporally bin the provided movie M by bin_factor.
%   M: Movie to be temporally binned [height x width x num_frames]
%   bin_factor: Number of frames in the original movie that will 
%       correspond to one frame in the binned movie.
%
% 2015 02 03 Tony Hyun Kim

[height, width, num_frames] = size(M);
num_downsampled_frames = floor(num_frames/bin_factor); % Note truncation
fprintf('Actual bin factor: %.3f\n', num_frames/num_downsampled_frames);

M_b = zeros(height, width, num_downsampled_frames, 'single');
for k = 1:num_downsampled_frames
    if (mod(k, 1000)==0)
        fprintf('  Frames %d / %d downsampled\n', k, num_downsampled_frames);
    end
    
    % Compute indices to original movie
    frame_start = 1 + bin_factor*(k-1);
    frame_end   = frame_start + (bin_factor-1);
    frames = frame_start:frame_end;
    
    % Compute the mean
    M_b(:,:,k) = mean(M(:,:,frames),3);
end
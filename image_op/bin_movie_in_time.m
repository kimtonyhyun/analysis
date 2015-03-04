function M_b = bin_movie_in_time(M, bin_factor, trial_indices)
% Temporally bin the provided movie M by bin_factor, but do not bin frames
% across trial boundaries.
%
% Inputs:
%   - M: Movie to be temporally binned [height x width x num_frames]
%   - bin_factor: Number of frames in the original movie that will 
%       correspond to one frame in the binned movie.
%   - trial_indices: [num_trials x 2] matrix indicating the start and end
%       frames of each trial with respect to the input movie M.
%

[height, width, num_frames] = size(M);

% Compute the number of binned frames for each trial
frames_per_trial = trial_indices(:,end) - trial_indices(:,1) + 1;
binned_frames_per_trial = floor(frames_per_trial/bin_factor);
num_binned_frames = sum(binned_frames_per_trial);

fprintf(['  bin_movie_in_time: %d frames from the original movie will be '...
         'omitted in binning to preserve trial boundaries.\n'],...
         num_frames - bin_factor*num_binned_frames);

M_b = 0;
% M_b = zeros(height, width, num_binned_frames, 'single');

% for k = 1:num_downsampled_frames
%     if (mod(k, 1000)==0)
%         fprintf('  Frames %d / %d downsampled\n', k, num_downsampled_frames);
%     end
%     
%     % Compute indices to original movie
%     frame_start = 1 + bin_factor*(k-1);
%     frame_end   = frame_start + (bin_factor-1);
%     frames = frame_start:frame_end;
%     
%     % Compute the mean
%     M_b(:,:,k) = mean(M(:,:,frames),3);
% end
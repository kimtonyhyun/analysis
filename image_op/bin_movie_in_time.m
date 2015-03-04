function [M_b, binned_indices] = bin_movie_in_time(M, bin_factor, trial_indices)
% Temporally bin the provided movie M by bin_factor, but do not bin frames
% across trial boundaries. Dangling frames will be excised. Furthermore,
% non-trial frames -- if they exist in the original movie -- will also be
% excised.
%
% Inputs:
%   - M: Movie to be temporally binned [height x width x num_frames]
%   - bin_factor: Number of frames in the original movie that will 
%       correspond to one frame in the binned movie.
%   - trial_indices: [num_trials x 4] matrix, where the i-th row indicates
%       the frame indices of Trial i as [start open-gate close-gate end]
%
% Outputs:
%   - M_b: Binned movie
%   - binned_indices: [num_trials x 4] matrix where the i-th row indicates
%       the frame indices with respect to the binned movie.
%

[height, width, ~] = size(M);

% Get the binned indices
binned_indices = bin_frame_indices(trial_indices, bin_factor);
binned_frames_per_trial = binned_indices(:,end) - binned_indices(:,1) + 1;
num_binned_frames = sum(binned_frames_per_trial);

% Generate the binned movie
M_b = zeros(height, width, num_binned_frames, 'single');

write_idx = 1;
for trial_idx = 1:num_trials
    for k = 1:binned_frames_per_trial(trial_idx)
        % Indices into the original movie
        frame_start = trial_indices(trial_idx,1) + bin_factor*(k-1);
        frame_end   = frame_start + (bin_factor-1);
        frames = frame_start:frame_end;
        
        % Compute the mean
        M_b(:,:,write_idx) = mean(M(:,:,frames),3);
        write_idx = write_idx + 1;
    end
end

end % bin_movie_in_time

function binned_indices = bin_frame_indices(orig_indices, bin_factor)
% TODO: Replace the globally visible 'io/bin_frame_indices.m' with this
%   implementation.

[num_trials, K] = size(orig_indices);

frames_per_trial = orig_indices(:,end) - orig_indices(:,1) + 1;
binned_frames_per_trial = floor(frames_per_trial/bin_factor);

% Explicitly remove dangling frames from the original frame indices
for trial_idx = 1:num_trials
    max_index = orig_indices(trial_idx,1) + ...
                bin_factor * binned_frames_per_trial(trial_idx) - 1;
    orig_indices(trial_idx,:) = min(orig_indices(trial_idx,:),...
                                    max_index*ones(1,K)); % Apply clamp
end

offsets = diff(orig_indices,1,2);
offsets = cumsum(offsets,2);

% Generate the binned indices
binned_start_indices = cumsum([1; binned_frames_per_trial(1:end-1)]);
binned_remaining_indices = repmat(binned_start_indices, 1, K-1) +...
                           floor(offsets/bin_factor);

binned_indices = [binned_start_indices binned_remaining_indices];

end
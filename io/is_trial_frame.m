function val = is_trial_frame(trial_frame_indices, test_frame_index)
% Tests if 'test_frame_index' is a trial frame
%
% Inputs:
%   trial_frame_indices: [num_trials x 2] matrix that indicates the
%       initial and final frame counts of each trial. The trials are
%       expected to be sorted in increasing order of the starting
%       frame number
%   test_frame_index: A frame number to test
%
% Example usage:
%   is_trial_frame([1 10; 15 20], 8): will return true

trial_idx = find(test_frame_index >= trial_frame_indices(:,1), 1, 'last');
val = test_frame_index <= trial_frame_indices(trial_idx,2);
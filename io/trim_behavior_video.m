function trim_behavior_video(plusmaze_source, behavior_source, trim)
% Extract the trial frames of the behavior video, so that the resulting
%   video lines up with the concatenated Miniscope movie (e.g. as produced
%   by `concatenateHDF5` or `concatenate_bigtiff`). The `trim` parameter
%   needs to match the values used in the Miniscope concatenation.
%
% Inputs:
%   plusmaze_source: Text file output from the plus maze
%   behavior_source: Behavior video (MPEG-4)
%   trim: Number of frames to drop from the beginning and end of each trial
%
% Example usage:
%   trim_behavior_video('mouse7_d07_ego-left.txt', 'mouse7_day07_ego-left.m4v', [15 5])
%
% Note: When the behavior video has fewer frames than the PlusMaze text
%   file, the behavior video will be uniformly "expanded" to match the
%   number of frames as indicated by the PlusMaze output.


% Get the trial frames (beginning and end) according to plusmaze output and
% trim the beginning and end
frame_indices = get_trial_frame_indices(plusmaze_source);
num_trials = size(frame_indices,1);
num_trial_frames = frame_indices(num_trials,4); % Very last frame
fprintf('  PlusMaze output (%s) has %d frames\n',...
    plusmaze_source, num_trial_frames);

% Frames to keep
trial_frame_indices = [frame_indices(:,1)+trim(1) frame_indices(:,4)-trim(2)];

% Expected number of trimmed frames
num_trimmed_frames = sum(diff(trial_frame_indices,1,2)+1);

% Read behavior video
behavior_video = VideoReader(behavior_source);
num_behavior_frames = behavior_video.NumberOfFrames;
fprintf('  Behavior video (%s) has %d frames\n',...
    behavior_source, num_behavior_frames);

if (num_trial_frames < num_behavior_frames)
    fprintf('  Behavior video has %d more frames than PlusMaze output. Aborting!\n',...
        num_behavior_frames - num_trial_frames);
else
    fprintf('  Behavior video is missing %d frames\n',...
        num_trial_frames - num_behavior_frames);
    frame_factor = double(num_behavior_frames) / double(num_trial_frames);
    fprintf('  Frame conversion factor of %.4f will be applied!\n',...
        frame_factor);

    % Prepare output file
    [~, name] = fileparts(behavior_source);
    output_name = sprintf('%s_trim%d-%d', name, trim(1), trim(2));
    trimmed_behavior_video = VideoWriter(output_name, 'MPEG-4');
    trimmed_behavior_video.FrameRate = 20; % FIXME: Don't hardcode
    
    write_idx = 0; % For tracking progress
    open(trimmed_behavior_video);
    for trial_idx = 1:num_trials
        for frame_idx = trial_frame_indices(trial_idx,1):trial_frame_indices(trial_idx,2)
            behavior_frame_idx = floor(frame_factor * frame_idx);
            A = read(behavior_video, behavior_frame_idx);
            A = rgb2gray(A);
            writeVideo(trimmed_behavior_video, A);
            
            write_idx = write_idx + 1;
            if (mod(write_idx,1000)==0)
                fprintf('%s: Frames %d of %d written\n',...
                    datestr(now), write_idx, num_trimmed_frames);
            end
        end
    end
    close(trimmed_behavior_video);
end
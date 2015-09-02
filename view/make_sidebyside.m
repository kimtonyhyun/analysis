function make_sidebyside(maze_source, behavior_source, miniscope_source, trials_to_record)

% Get the table of frame indices
%------------------------------------------------------------
frame_indices = get_trial_frame_indices(maze_source);
frame_indices = frame_indices(:,[1 end]); % [Trial-start Trial-end]
num_frames = frame_indices(end,end);

% Open the behavior video and make sure the number of frames matches the
% trial index table
%------------------------------------------------------------
behavior_vid = VideoReader(behavior_source);
behavior_width = behavior_vid.Width;
behavior_height = behavior_vid.Height;
assert(num_frames == behavior_vid.NumberOfFrames,...
    sprintf('Number of frames in behavior video (%d) does not match that in the maze file (%d)!\n',...
            behavior_vid.NumberOfFrames, num_frames));

% Open the miniscope movie and make sure that the number of frames matches
% the trial index table
%------------------------------------------------------------
movie_size = get_movie_info(miniscope_source); % [h, w, num_frames]
assert(num_frames == movie_size(3),...
    sprintf('Number of frames in miniscope recording (%d) does not match that in the maze file (%d)!\n',...
            movie_size(3), num_frames));
        
% Open one trial of the miniscope movie, to compute the appropriate CLim
M = load_movie_from_hdf5(miniscope_source, frame_indices(1,:));
movie_clim = compute_movie_scale(M);
movie_clim_range = diff(movie_clim);

% Also, determine the width when the miniscope frame is rescaled to match
% the behavior frame height
m = imresize(M(:,:,1), [behavior_height NaN]);
A = zeros(behavior_height, behavior_width + size(m,2), 'single');
h = imagesc(A, [0 1]); % Images will be scaled to [0 1]
colormap gray;
truesize;

% Prepare the output file
%------------------------------------------------------------
[~, behavior_base_name, ~] = fileparts(behavior_source);
output_name = sprintf('%s_sbs.avi', behavior_base_name);
fprintf('%s: Side-by-side video will be saved to "%s"...\n', datestr(now), output_name);

writerObj = VideoWriter(output_name, 'Grayscale AVI');
writerObj.FrameRate = 20; % FIXME
open(writerObj);

for trial_idx = trials_to_record
    fprintf('%s: Processing trial %d...\n', datestr(now), trial_idx);
    
    trial_frames = frame_indices(trial_idx,:);
    num_frames_in_trial = diff(trial_frames) + 1;
    M_trial = load_movie_from_hdf5(miniscope_source, frame_indices(trial_idx,:));
    
    for k = 1:num_frames_in_trial
        b = rgb2gray(behavior_vid.read(trial_frames(1) + (k-1))); % uint8
        b = single(b) / 255; % Rescaled to [0 1]
        
        m = imresize(M_trial(:,:,k), [behavior_height NaN]); % Match to behavior height
        m = (m-movie_clim(1)) / movie_clim_range; % Rescaled to [0 1]
        
        sbs = [b m];
        sbs = max(sbs, 0); % Clamp at 0 from below
        sbs = min(sbs, 1); % Clamp at 1 from top
        
        set(h, 'CData', sbs);
        drawnow;
        
        writeVideo(writerObj, sbs);
    end
end
close(writerObj);
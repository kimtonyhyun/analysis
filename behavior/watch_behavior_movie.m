function watch_behavior_movie( movie, trial, plusmaze_txt )
% Plays behavior video for specified trial
%
% Example use:
% trials = [12:20 24 26];
% watch_behavior_movie('c9m7d09.mp4',trials,'c9m7d09.txt')
% 
% Input:
%     movie: behavior video
%     trial: matrix of trial(s) to watch, e.g. [1 4 5 8]
%     plusmaze_txt: plusmaze text file, or trimmed/binned version
%
% 2015-03-13 Fori Wang
%
%

    % parse plusmaze txt file
    trial_indices = plusmaze_txt;
    trial_indices = get_trial_frame_indices(trial_indices);
    trial_indices = [trial_indices(:,1) trial_indices(:,4)];

    % read in behavior movie
    fprintf('Reading in Behavior Movie...\n');
    behavior_vid = VideoReader(movie);
    
    % initialize plot
    figure;
    
    this_trial_indices = trial_indices(1,:);
    trial_start = this_trial_indices(1);  
    first_image = read(behavior_vid,trial_start);
    k = imagesc(first_image);
    title(sprintf('Trial %d',trial(1)),'FontSize',18);
    axis image; colormap gray;
    
    blank_image = zeros(size(first_image));
    
    for i = 1:length(trial)
        
        % update trial indices
        this_trial = trial(i);
        title(sprintf('Trial %d',this_trial));
        this_trial_indices = trial_indices(this_trial,:);
        trial_start = this_trial_indices(1);  
        
        % load trial video
        trial_video = read(behavior_vid,this_trial_indices);
        trial_video = squeeze(trial_video(:,:,1,:)); % 3D movie stack
       
        % play trial frames
        for frame_idx = 1:size(trial_video,3);
            image = trial_video(:,:,frame_idx);
            set(k,'CData',image);          
            pause(0.01);
        end
        
        % display black screen at end of each trial
        set(k,'CData',blank_image);
        pause(0.1);
    end

end

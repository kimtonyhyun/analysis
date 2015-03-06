function watch_trial_with_centroids( movie, centroids_mat, trial, varargin )
% Plays behavior video with centroid on top for specified trial
% 
% Input:
%     movie: behavior video
%     centroids_mat: position data from get_mouse_XY_pos
%     trial: matrix of trial(s) to watch, e.g. [1 4 5 8]
%     varargin: plusmaze txt file, e.g. 'mouse7_day09_allo-south.txt' or
%               trimmed/binned version
%
% Example use:
% trials = 1:15;
% watch_trial_with_centroids('mouse7.mp4','centroids.mat',trials,'mouse7.txt')
%
% 2015-03-06 Fori Wang
%
%

    % get trial indices
    if ~isempty(varargin)
        trial_indices = varargin{1};
        trial_indices = get_trial_frame_indices(trial_indices);
        trial_indices = [trial_indices(:,1) trial_indices(:,4)];
    else
        [trial_indices,~] = find_start_end_of_trials(centroids_mat);
    end

    % read in behavior movie
    fprintf('Reading in Behavior Movie...\n');
    behavior_vid = VideoReader(movie);
    
    % load centroids
    load(centroids_mat)
    
    % initialize plot
    figure;
    
    this_trial_indices = trial_indices(1,:);
    trial_start = this_trial_indices(1);  
    first_image = read(behavior_vid,trial_start);
    k = imagesc(first_image);
    title('Trial');
    axis image; colormap gray;
    hold on
    l = plot(centroids(trial_start,1),centroids(trial_start,2),'g*');
    
    blank_image = zeros(size(first_image));
    
    for i = 1:length(trial)
        
        % update trial indices
        title(sprintf('Trial %d',trial(i)));
        this_trial_indices = trial_indices(i,:);
        trial_start = this_trial_indices(1);  
        
        % load trial video
        trial_video = read(behavior_vid,this_trial_indices);
        trial_video = squeeze(trial_video(:,:,1,:)); % 3D movie stack

        % update plot
        for frame_idx = 1:size(trial_video,3);
            image = trial_video(:,:,frame_idx);
            this_centroid = centroids(frame_idx+trial_start-1,:);

            set(k,'CData',image);
            set(l,'XData',this_centroid(1),'YData',this_centroid(2));
            pause(0.1);
        end
        
        % display black screen at end of each trial
        set(k,'CData',blank_image);
        pause(0.5);
    end

end


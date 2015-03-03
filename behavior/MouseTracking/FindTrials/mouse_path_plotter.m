function mouse_path_plotter( movie, centroids_mat)
% Plots the path of the mouse during each trial based on position data and
% trial frames extracted from behavior movie
% 
% Input:
%     movie: Behavior video
%     centroids_mat: position data from get_mouse_XY_pos
%
% 
% 2015-03-01 Fori Wang

    % find trial start and ends
    [trial_frames,~]=find_start_end_of_trials(centroids_mat);
    
    % load centroids
    load(centroids_mat)% FIX: assumes matrix saved as 'centroids'
    
    % read in behavior movie
    fprintf('Reading in Behavior Movie...');
    behavior_vid = VideoReader(movie);
    
    figure;
    
    for trial_idx = 1:length(trial_frames)
        trial_start = trial_frames(trial_idx,1);
        trial_end = trial_frames(trial_idx,2);
        scrollsubplot(5,6,trial_idx)
        
        %set image as first image of trial
        image = read(behavior_vid,trial_start);
        imagesc(image);
        title(sprintf('Trial %d',trial_idx));
        axis image;
        set(gca,'XTick',[],'YTick',[]);
        colormap gray;
        hold on
        plot(centroids(trial_start:trial_end,1),...
             centroids(trial_start:trial_end,2),...
             'g-');
    end

end
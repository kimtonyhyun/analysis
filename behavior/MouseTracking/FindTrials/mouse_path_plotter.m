function mouse_path_plotter( movie, centroids_mat, trial_frames )
% Plots the path of the mouse during each trial based on position data and
% trial frames extracted from behavior movie
% 
% Input:
%     movie: Behavior video
%     centroids_mat: position data from get_mouse_XY_pos
%     trial_frames: trial_frames from find_start_end_of_trials
% 
% 2015-03-01 Fori Wang
    
    load(centroids_mat)% FIX: assumes matrix saved as 'centroids'
    figure;
    behavior_vid = VideoReader(movie);
    image = read(behavior_vid,[1 5000]);
    image = squeeze(image(:,:,1,:));
    avg_image = mean(image,3);
%     imagesc(avg_image);
%     title('Mouse Path Plotter')
%     axis image
%     hold on
    
%     color_options = hsv(length(trial_frames));
    
    for trial_idx = 31:60 %length(trial_frames)
        trial_start = trial_frames(trial_idx,1);
        trial_end = trial_frames(trial_idx,2);
        subplot(5,6,trial_idx-30)
        imagesc(avg_image);
        axis image
        hold on
        plot(centroids(trial_start:trial_end,1),...
             centroids(trial_start:trial_end,2),...
             'g-');
%              '-','Color',color_options(trial_idx,:));
    end

end


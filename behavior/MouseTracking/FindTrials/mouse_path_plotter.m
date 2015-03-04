function [trial_frames] = mouse_path_plotter( movie, centroids_mat, varargin)
% Plots the path of the mouse during each trial based on position data and
% trial frames (either specified by plusmaze txt or extracted from movie
% using find_start_end_of_trials)
% 
% Input:
%     movie: Behavior video
%     centroids_mat: position data from get_mouse_XY_pos
%     varargin: plusmaze txt file, e.g. 'mouse7_day09_allo-south.txt'
% 
% Output:
%     trial_frames: if no plusmaze source, can return trial_frames generated
%         by find_start_end_of_trials
% 
% 2015-03-01 Fori Wang

    % read in plusmaze source if available
    if ~isempty(varargin)
        maze_data = varargin{1};
        [frame_indices,~,~]=parse_plusmaze(maze_data);
        trial_frames = [frame_indices(:,1) frame_indices(:,4)];
    else
        % if no plusmaze source, find trial start and ends from movie
        [trial_frames,~]=find_start_end_of_trials(centroids_mat);
    end
    
    % load centroids
    load(centroids_mat)% assumes matrix saved as 'centroids'
    
    % read in behavior movie
    fprintf('Reading in Behavior Movie...');
    behavior_vid = VideoReader(movie);
    
    figure;
    
    for trial_idx = 1:size(trial_frames,1)
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
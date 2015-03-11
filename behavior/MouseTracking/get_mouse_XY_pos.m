function [ centroids ] = get_mouse_XY_pos( movie, varargin )
% Extract the X,Y positions ( centroids ) of the mouse from the behavior
%   video ( movie ), with live monitor option (specify 'displayTracking')
%  
% Example uses:
%     Normal mode:
%       centroids = get_mouse_XY_pos('c9m7d08_ti2-sub.mp4');
%   Live monitor mode:	
%         centroids = get_mouse_XY_pos('c9m7d08_ti2-sub.mp4','displayTracking');
%     Trial-aware mode:
%       centroids = get_mouse_XY_pos('c9m7d08_ti2-sub.mp4','c9m7d08_cr_ti2.txt');
% 
% Input:
%   - movie: Behavior video (m4v); movie should be cropped (wavy curtain
%   leads to fall centroids (haven't found another way to fix this)
%   - varargin:
%     'displayTracking': Displays side-by-side of original movie and
%           processing movie, along with centroid on each. blue * means good
%           centroid value, red * means could not find mouse, using previous
%           centroid value, magenta * means temporary centroid assigned at
%           trial boundary (to be reassigned same as next centroid)
%     'trial_indices.txt': If provided, will not use prev_centroid method
%     at the trial boundaries to avoid bleed-through across trials
%
% Returns matrix centroids where each row is an (x,y) coordinate of the
%   mouse. 1st column = X, 2nd column = Y. Length of matrix =
%   number of frames in movie
%
% Notes: Analyzes movie in chunks of chunk_size=1000 frames (smaller memory load)
%
% Updated 2015-03-10 Fori Wang

    display_tracking = 0;
    % check if live monitor option requested
    if ~isempty(varargin)
        option = varargin{1};
        if strcmp(option,'displayTracking')
            display_tracking = 1;
            figure;
            fprintf('Displaying live tracking...');            
        else
            trial_aware = 1;
            trial_indices = get_trial_frame_indices(option);
            trial_indices = trial_indices(:,[1 4]);
        end
    end
        
    behavior_vid = VideoReader(movie);
    num_frames = behavior_vid.NumberofFrames;

    % initialize centroids
    centroids = zeros(num_frames,2);
    
    % initialize variables for trial_aware mode
    lost_centroids=0;
    reassign_prev_centroid=0;
    
    % setup chunks of frames to read in movie
    chunk_size = 1000;
    frame_chunks = make_frame_chunks(num_frames,chunk_size);
    
    % setup background image: average of 5000 frames
    bg_vid = read(behavior_vid,[1 5000]);
%     bg_vid = squeeze(bg_vid(:,:,1,:)); % 3D movie stack
    bg_image = mean(bg_vid,4);
    bg_image = uint8(bg_image);
    
    c_old = [0 0];
    
    for idx = 1:length(frame_chunks)

        fprintf('  Processing frames %d to %d (%s)\n',...
            frame_chunks(idx,1), frame_chunks(idx,2), datestr(clock));

        % read in all of the frames for the trial at the beginning
        frame_range = frame_chunks(idx,:);
        video = read(behavior_vid,frame_range);
%         video = squeeze(video(:,:,1,:)); % 3D movie stack
        
        if display_tracking
            
            % plot original image on the left
            subplot(121);
            image = video(:,:,:,1);
            h = imagesc(image);
            title(sprintf('Original: Frames %d - %d',...
                          frame_chunks(idx,1), frame_chunks(idx,2)));
            axis image; colormap gray;
            hold on
            
            % plot initial centroids on original video
            k = plot(0,0,'k*');

            % plot initial tracking image on the right
            subplot(122);
            j = imagesc(image); colormap gray;
            title(sprintf('Tracker: Frames %d - %d',...
                          frame_chunks(idx,1),frame_chunks(idx,2)));
            hold on
            axis image;

            % plot centroids on tracking video
            l = plot(0,0,'k*');
        end

        for frame_idx = 1:size(video,4)
           
           % Update original image CData
            image = video(:,:,:,frame_idx);
%             image = im2double(image,'indexed');
            if display_tracking
                set(h,'CData',image);
                pause(0.00001);
            end
           
            % Find the mouse blob using findMouse helper function
            thresh = 20; % assumes that black mouse has RGB values <20
            [final_image,s] = findMouse(bg_image,image,thresh);
            area_vector = [s.Area];
            length_vector = [s.MinorAxisLength];
            
            if display_tracking
                % Update tracking
                set(j,'CData',final_image);
                pause(0.00001);
            end

            % Save centroids data and update plot
            true_frame_idx = frame_chunks(idx,1)+frame_idx-1;           
                
            if isempty(area_vector) %sometimes blob disappears, go to previous centroid
                
                % trial_aware mode
                % if image is the first image of a trial, use next centroid
                if trial_aware && find(true_frame_idx==trial_indices(:,1),1)
                    reassign_prev_centroid = 1;
                    lost_centroids = lost_centroids+1;
                    centroid_color = 'm';
                elseif reassign_prev_centroid
                    % in cases mouse also lost in frames immediately
                    % following first frame
                    lost_centroids = lost_centroids+1;
                    centroid_color = 'm';
                else % use previous centroid
                    centroids(true_frame_idx,:) = c_old;
                    centroid_color = 'r';
                end
                
            else % new centroid
                [~, id] = max(length_vector); %assume mouse is the fattest blob
                c_new = s(id(1)).Centroid(1:2);
                centroids(true_frame_idx,:) = c_new;
                
                % trial_aware mode
                if reassign_prev_centroid % assign c_new to previous centroid where blob was lost
                    centroids(true_frame_idx-lost_centroids:true_frame_idx-1,:) = c_new;
                    reassign_prev_centroid = 0;
                    lost_centroids = 0;
                end
                
                c_old = c_new;
                centroid_color = 'b';
            end
            
            if display_tracking
                % Update subplots
                set(k,'XData',c_old(1),'YData',c_old(2),'color',centroid_color);
                set(l,'XData',c_old(1),'YData',c_old(2),'color',centroid_color);
                pause(0.00001);
            end
        end    
    end
    
    centroids(1,:)=centroids(2,:); %get rid of c_old [0 0] start

end


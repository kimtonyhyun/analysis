function [ centroids ] = get_mouse_XY_pos( movie, varargin )
% Extract the X,Y positions ( centroids ) of the mouse from the behavior
%   video ( movie ), with live monitor option (specify 'displayTracking')
%  
% Input:
%   - movie: Behavior video (m4v); movie should be cropped (wavy curtain
%   leads to fall centroids (haven't found another way to fix this)
%   - varargin: 'displayTracking' displays side-by-side of original movie and
%   processing movie, along with centroid on each. blue * means good
%   centroid value, red * means could not find mouse, using previous
%   centroid value
%
% Returns matrix centroids where each row is an (x,y) coordinate of the
%   mouse. 1st column = X, 2nd column = Y. Length of matrix =
%   number of frames in movie
%
% Notes: Analyzes movie in chunks of chunk_size=1000 frames (smaller memory load)
%   *** read function does not work properly on Linux, so unfortunately
%   have to run this on a Windows machine ***
%   *** Not very fast, but only have to run once; ~54K frames = ~40 min
%
% 2015-02-26 Fori Wang

    display_tracking = 0;
    % check if live monitor option requested
    if ~isempty(varargin)
        option = varargin{1};
        
        if isfield(options,'displayTracking')
            display_tracking=1;
        else
            fprintf('Cannot recognize this option. Did you mean displayTracking?')
        end
    end
        
    behavior_vid = VideoReader(movie);
    num_frames = behavior_vid.NumberofFrames;

    % initialize centroids
    centroids = zeros(num_frames,2);
%     centroids = cell(num_frames,1);
    
    % setup chunks of frames to read in movie
    chunk_size = 1000;
    frame_chunks = make_frame_chunks(num_frames,chunk_size);
    
    % setup background image: average of 5000 frames
    bg_vid = read(behavior_vid,[1 5000]);
    bg_vid = squeeze(bg_vid(:,:,1,:)); % 3D movie stack
    bg_image = mean(bg_vid,3);
%     bg_image = uint8(bg_image);
    
    c_old = [0 0];
    
    for idx = 1:length(frame_chunks)

        fprintf('  Processing frames %d to %d (%s)\n',...
            frame_chunks(idx,1), frame_chunks(idx,2), datestr(clock));

        % read in all of the frames for the trial at the beginning
        frame_range = frame_chunks(idx,:);
        video = read(behavior_vid,frame_range);
        video = squeeze(video(:,:,1,:)); % 3D movie stack
        
        if display_tracking
            
            % plot original image on the left
            subplot(121);
            image = video(:,:,:,1);
            h = imagesc(image);
            title(sprintf('Original: Frames %d - %d',...
                          frame_indices(idx,1), frame_indices(idx,2)));
            axis image; 
            hold on
            
            % plot initial centroids on original video
            k = plot(0,0,'k*');

            % plot initial tracking image on the right
            subplot(122);
            j = imagesc(image);
            title(sprintf('Tracker: Frames %d - %d',...
                          frame_indices(idx,1),frame_indices(idx,2)));
            hold on
            axis image;

            % plot centroids on tracking video
            l = plot(0,0,'k*');
        end

        for frame_idx = 1:size(video,3)
           
           % Update original image CData
            image = video(:,:,frame_idx);
            image = im2double(image,'indexed');
            if display_tracking
                set(h,'CData',image); 
            end
           
            % Find the mouse blob using findMouse helper function
            thresh = 20; % assumes that black mouse has RGB values <20
            [final_image,s] = findMouse(bg_image,image,thresh);
            area_vector = [s.Area];
            length_vector = [s.MinorAxisLength];
            
            if display_tracking
                % Update tracking
                set(j,'CData',final_image);
            end

            % Save centroids data and update plot
            if isempty(area_vector) %sometimes blob disappears, go to previous centroid
                centroids(frame_chunks(idx,1)+frame_idx-1,:) = c_old;
                
                if display_tracking
                    % Update subplots
                    set(k,'XData',c_old(1),'YData',c_old(2),'color','r');
                    set(l,'XData',c_old(1),'YData',c_old(2),'color','r');
                end
                
            else % new centroid
                [~, id] = max(length_vector); %assume mouse is the fattest blob
                c_new = s(id(1)).Centroid(1:2);
                centroids(frame_chunks(idx,1)+frame_idx-1,:) = c_new;
                c_old = c_new;
                
                if display_tracking
                    % Update subplots
                    set(k,'XData',c_new(1),'YData',c_new(2),'color','b');
                    set(l,'XData',c_new(1),'YData',c_new(2),'color','b');
                end
            end
        end    
    end
    
    centroids(1,:)=centroids(2,:); %get rid of c_old [0 0] start

end


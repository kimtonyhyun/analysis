function [ centroids ] = get_mouse_XY_pos_live_monitor( movie )
% Extract the X,Y positions ( centroids ) of the mouse from the behavior
%   video ( movie )
%   
% Input:
%   movie: Behavior video (m4v); movie should be cropped (wavy curtain
%   leads to fall centroids (haven't found another way to fix this)
%
% Returns matrix centroids where each row is an (x,y) coordinate of the
%   mouse. 1st column = X, 2nd column = Y. Length of matrix =
%   number of frames in movie
%
% Notes: Analyzes movie in chunks of 500 frames (smaller memory load)
%   *** read function does not work properly on Linux, so unfortunately
%   have to run this on a Windows machine ***
%   *** Not very fast, but only have to run once; ~54K frames = ~40 min
%
% 2015-02-26 Fori Wang

    behavior_vid = VideoReader(movie);
    num_frames = behavior_vid.NumberofFrames;

    % initialize centroids
    centroids = cell(num_frames,1);
    
    % setup chunks of frames to read in movie
    frame_indices = make_frame_indices(num_frames,1000);
    
    % setup background image: average of 5000 frames
    bg_vid = read(behavior_vid,[1 5000]);
    bg_image = mean(bg_vid,4);
    bg_image = uint8(bg_image);
    
    for idx = 1:length(frame_indices)

        fprintf('  Processing frames %d to %d (%s)\n',...
            frame_indices(idx,1), frame_indices(idx,2), datestr(clock));

        % read in all of the frames for the trial at the beginning
        frame_range = [frame_indices(idx,1) frame_indices(idx,2)];
        video = read(behavior_vid,frame_range);

        
              
        % plot original image (original movie on the left)
        subplot(121);
        image = video(:,:,:,1);
        h = imagesc(image);
        title(sprintf('Original: Frames %d - %d',...
                      frame_indices(idx,1), frame_indices(idx,2)));
        axis image; 
        hold on
        c_old = [0 0];

        % plot initial centroids on Original Video
        k = plot(0,0,'k*');

        % plot initial Tracker image (on the right)
        subplot(122);
        j = imagesc(image);
        title(sprintf('Tracker: Frames %d - %d',...
                      frame_indices(idx,1),frame_indices(idx,2)));
        hold on
        axis image;

        % plot centroids on Tracker subplot
        l = plot(0,0,'k*');
        

        for frame_idx = frame_indices(idx,1)-(frame_indices(idx,1)-1):...
                        frame_indices(idx,2)-(frame_indices(idx,1)-1)
           
           % Update original image CData
            image = video(:,:,:,frame_idx);
            set(h,'CData',image); 
           
            % Find the mouse blob using findMouse helper function
            thresh = 20; % assumes that black mouse has RGB values <20
            [final_image,s] = findMouse(bg_image,image,thresh);
            area_vector = [s.Area];
            length_vector = [s.MinorAxisLength];
            
            % Update binarized image
            set(j,'CData',final_image);

            % Save centroids data and update plot
            if isempty(area_vector) %sometimes blob disappears, go to previous centroid
                centroids{frame_indices(idx,1)+frame_idx-1} = c_old;

                % Update subplots
                set(k,'XData',c_old(1),'YData',c_old(2),'color','r');
                set(l,'XData',c_old(1),'YData',c_old(2),'color','r');
                
            else % new centroid
                [~, id] = max(length_vector); %assume mouse is the fattest blob
                c_new = s(id(1)).Centroid(1:2);
                centroids{frame_indices(idx,1)+frame_idx-1} = c_new;
                c_old = c_new;

                % Update subplots
                set(k,'XData',c_new(1),'YData',c_new(2),'color','b');
                set(l,'XData',c_new(1),'YData',c_new(2),'color','b');
                
            end
        end    
    end
    
    centroids=cell2mat(centroids); %convert to matrix
    centroids(1,:)=centroids(2,:); %get rid of c_old [0 0] start
    
end


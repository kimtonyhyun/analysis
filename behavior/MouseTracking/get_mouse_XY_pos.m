function [ centroids ] = get_mouse_XY_pos( movie )
% Extract the X,Y positions ( centroids ) of the mouse from the behavior
%   video ( movie )
%   
% Input:
%   movie: Behavior video (m4v)
%
% Returns cell (centroids) where each row is an (x,y) coordinate of the
%   mouse. Length of cell = number of frames in moive
%
% Notes: Analyzes movie in chunks of 500 frames (smaller memory load)
%
% 2015-02-26 Fori Wang

    behavior_vid = VideoReader(movie);
    num_frames = behavior_vid.NumberofFrames;

    % initialize centroids
    centroids = cell(num_frames,1);

    frame_indices = make_frame_indices(num_frames);

    for idx = 1:length(frame_indices)

        fprintf('  Processing frames %d to %d (%s)\n',...
            frame_indices(idx,1), frame_indices(idx,2), datestr(clock));

        % read in all of the frames for the trial at the beginning
        frame_range = [frame_indices(idx,1) frame_indices(idx,2)];
        video = read(behavior_vid,frame_range);       
        
%         % plot original image (original movie on the left)
%         subplot(121);
%         image = video(:,:,:,1);
%         h = imagesc(image);
%         title(sprintf('Original: Frames %d - %d',...
%                       frame_indices(idx,1), frame_indices(idx,2)));
%         axis image; 
%         hold on
        c_old = [0 0];

%         % plot initial centroids on Original Video
%         k = plot(0,0,'k*');

%         % plot initial Tracker image (on the right)
%         subplot(122);
%         j = imagesc(image);
%         title(sprintf('Tracker: Frames %d - %d',...
%                       frame_indices(idx,1),frame_indices(idx,2)));
%         hold on
%         axis image;

%         % plot centroids on Tracker subplot
%         l = plot(0,0,'k*');

        % set background subtraction image to be first frame
        bg_image = video(:,:,:,frame_indices(idx,2)-(frame_indices(idx,1)-1));

        for frame_idx = frame_indices(idx,1)-(frame_indices(idx,1)-1):...
                        frame_indices(idx,2)-(frame_indices(idx,1)-1)
           fprintf('  Processing frame %d until %d (%s)\n',...
                frame_indices(idx,1)+frame_idx-1,...
                frame_indices(idx,2), datestr(clock));

           % Update original image CData
            image = video(:,:,:,frame_idx);
%             set(h,'CData',image); 
%             pause(0.001);

            % Find the mouse blob using findMouse helper function
            thresh = 20; % assumes that black mouse has RGB values <20
            [~,s] = findMouse(bg_image,image,thresh);
            area_vector = [s.Area];
            length_vector = [s.MinorAxisLength];

            % Update binarized image
%             set(j,'CData',final_image);
%             pause(0.001);

            % Save centroids data and update plot
            if isempty(area_vector) %sometimes blob disappears, go to previous centroid
                centroids{frame_indices(idx,1)+frame_idx-1} = c_old;

%                 % Update subplots
%                 set(k,'XData',c_old(1),'YData',c_old(2),'color','r');
%                 set(l,'XData',c_old(1),'YData',c_old(2),'color','r');
%                 pause(0.001);

            else % new centroid
                [~, id] = max(length_vector); %assume mouse is the fattest blob
                c_new = s(id(1)).Centroid(1:2);
                centroids{frame_indices(idx,1)+frame_idx-1} = c_new;
                c_old = c_new;

%                 % Update subplots
%                 set(k,'XData',c_new(1),'YData',c_new(2),'color','b');
%                 set(l,'XData',c_new(1),'YData',c_new(2),'color','b');
%                 pause(0.001);
            end
        end    
    end

end


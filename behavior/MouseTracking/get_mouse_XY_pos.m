function [ centroids ] = get_mouse_XY_pos( movie )
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

    frame_indices = make_frame_indices(num_frames);

    for idx = 1:length(frame_indices)

        fprintf('  Processing frames %d to %d (%s)\n',...
            frame_indices(idx,1), frame_indices(idx,2), datestr(clock));

        % read in all of the frames for the trial at the beginning
        frame_range = [frame_indices(idx,1) frame_indices(idx,2)];
        video = read(behavior_vid,frame_range);       
        
        c_old = [0 0];

        % set background subtraction image to be later time point % not a
        % great solution, but can't think of anything better
        bg_offset = 240; % ~12 seconds; for 5/5 sessions, 140 prob better
        if idx+240 > frame_indices(idx,2)
            bg_image = video(:,:,:,frame_indices(idx,1)-bg_offset);
        else
            bg_image = video(:,:,:,frame_indices(idx,1)+bg_offset);
%         bg_image = video(:,:,:,frame_indices(idx,2)-(frame_indices(idx,1)-1));

        for frame_idx = frame_indices(idx,1)-(frame_indices(idx,1)-1):...
                        frame_indices(idx,2)-(frame_indices(idx,1)-1)
           
           % Update original image CData
            image = video(:,:,:,frame_idx);

            % Find the mouse blob using findMouse helper function
            thresh = 20; % assumes that black mouse has RGB values <20
            [~,s] = findMouse(bg_image,image,thresh);
            area_vector = [s.Area];
            length_vector = [s.MinorAxisLength];

            % Save centroids data and update plot
            if isempty(area_vector) %sometimes blob disappears, go to previous centroid
                centroids{frame_indices(idx,1)+frame_idx-1} = c_old;

            else % new centroid
                [~, id] = max(length_vector); %assume mouse is the fattest blob
                c_new = s(id(1)).Centroid(1:2);
                centroids{frame_indices(idx,1)+frame_idx-1} = c_new;
                c_old = c_new;

            end
        end    
    end
    
    centroids=cell2mat(centroids) %convert to matrix

end


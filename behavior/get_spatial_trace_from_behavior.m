function spatial_trace = get_spatial_trace_from_behavior(sample_image,video)   
%function spatial_trace = get_spatial_trace_from_behavior(sample_image,video)
%
%Extract X and Y coordinates of the mouse from the behavior video
%
%   sample_image : A snapshot the behavior setup for use in the maze region
%   extraction step.
%
%   video : [Height x Width x num_frames] array
%
%   spatial_trace : [num_frames x 2] matrix where the 1st column is the X-coordinate and
%   the 2nd column is the Y-coordinate of the mouse. 

    
    num_frames = size(video,3);
    centroids = zeros(num_frames,2); 
    
    %anything above thresh_whiteness is considered white
    thresh_whiteness = 100;
    
    %anything below thresh_blackness is considered black
    thresh_blackness = 20;
    c_old = [0 0];
    
    %extract roughly the maze from the sample image (based on whiteness)
    if size(sample_image,3)>1
        maze_image = rgb2gray(sample_image); %convert to grayscale
    else
        maze_image = sample_image;
    end
    maze_image = medfilt2(maze_image,[20,20]); %get rid of details
    maze_image = maze_image > thresh_whiteness; %binarize with threshold
    for frame_idx = 1:num_frames
        % Update original image CData
        frame = video(:,:,frame_idx);  
        %make off-maze region bright (bright is non-interesting)
        frame(maze_image == 0) = 255;
        % binarize image by thresholding (black is interesting - mouse is black)
        binary_image = (frame < thresh_blackness) ;       
        %get rid of smaller dots that are not the mouse
        binary2 = bwareaopen(binary_image,100,8);      
        % get centroid, area, and length of blobs
        s = regionprops(logical(binary2(:,:,1)), 'area','centroid','minoraxislength');
        area_vector = [s.Area];
        length_vector = [s.MinorAxisLength];
        centroid_vector = [s.Centroid];
        %we need to reshape this to a num_centroids by 2 matrix:
        centroid_vector = reshape(centroid_vector,2,length(centroid_vector)/2)';
              
        if isempty(area_vector) %sometimes blob disappears, go to previous centroid
            centroids(frame_idx,:) = c_old;
        else %when the blob is there (should be most times)
            [~, idx] = max(length_vector); %assume mouse is the fattest blob
            c_new = centroid_vector(idx(1),:);
%             c_new = s(idx(1)).Centroid(1:2);
            centroids(frame_idx,:) = c_new;            
            % Update centroids on original and binary2 using new data
            c_old = c_new;
        end    
    end
    spatial_trace = centroids;

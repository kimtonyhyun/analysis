function centroid_plotter( movie, centroids_mat )
% Plots all centroids on first frame of behavior movie for general inspection
% 
% Input:
%     Behavior movie
%     centroids generated by get_mouse_XY_pos function
%     
% 2015-02-27 Fori Wang
    load(centroids_mat)% FIX: assumes matrix saved as 'centroids'
    figure;
    behavior_vid = VideoReader(movie);
    image = read(behavior_vid,1); %grab first frame
    imagesc(image)
    title('Centroid Plotter')
    axis image
    hold on
    
    plot(centroids(:,1),centroids(:,2),'g*');
end
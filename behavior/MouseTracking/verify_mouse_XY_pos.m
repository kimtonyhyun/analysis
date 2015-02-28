function verify_mouse_XY_pos( movie, centroids )
% Plots centroids obtained from get_mouse_XY_pos onto original behavior
%   video so that user can verify that centroids are accurate
%
% Input:
%   movie: Behavior video (m4v)
%
% Plots movie with * where centroid should be
%
% 2015-02-27 Fori Wang

    behavior_vid = VideoReader(movie);
    num_frames = behavior_vid.NumberofFrames;

    % plot initial image
    figure(1);
    image = read(behavior_vid,1);
    h = imagesc(image);
    title('MouseTracker');
    axis image;
    hold on
    
    % plot initial centroid
    j = plot(0,0,'g*');
        
    for frame_idx = 1:num_frames
        
        fprintf('  Processing frame %d of %d (%s)\n',...
                frame_idx,num_frames,datestr(clock));
                
        image = read(behavior_vid,frame_idx);
        set(h,'CData',image);
                
        this_centroid = centroids{frame_idx};        
        set(j,'XData',this_centroid(1),'YData',this_centroid(2)); 
        pause(0.0001);
    end
end


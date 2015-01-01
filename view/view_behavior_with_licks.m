function view_behavior_with_licks(video_source,lick_source)
%function view_behavior_with_licks(video_source,lick_source)
%
%View the behavior video with a superimposed rectangular box in the upper 
%left corner that indicates whether there is a lick at a given moment 
%
%   video_source: .mp4 file to be read as the behavior video source
%
%   lick_source: .txt file that contains a binary array with length equal 
%   to thenumber of frames in the video_source (syncing of the lick_source
%   and video_source is assumed)
%
%Hakan Inan(Dec 14)
%

% Read the lick file
fileID = fopen(lick_source,'r');
lick_data = fscanf(fileID,'%d');
fclose(fileID);
num_frames = length(lick_data);

% Read video file
behavior_vid = VideoReader(video_source);
h = imshow(read(behavior_vid,1));
video_source_name = strrep(video_source, '_', '\_');

% Display video with a block that indicates lick or no lick
block = [1 60 1 60]; % [x1 x2 y1 y2]
for k = 1:num_frames
    frame = read(behavior_vid,k);
    frame = frame(:,:,1);
    
    % Lick indicator
    frame(block(1):block(2),...
          block(3):block(4)) = 255; % Border
    if (~lick_data(k)) % No lick
        frame((block(1)+1):(block(2)-1),...
              (block(3)+1):(block(4)-1)) = 0; % Block interior is dark
    end
    
    % Update the figure
    set(h, 'CData', frame);
    title(sprintf('%s (Frame %d/%d)',...
                  video_source_name, k, num_frames));
    drawnow;
end



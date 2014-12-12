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

%Read the lick file
fileID = fopen(lick_source,'r');
lick_data = fscanf(fileID,'%d');
num_frames = length(lick_data);
% Read video file 
behavior_vid = VideoReader(video_source);
figure;
%the shape to display to indicate lick/no lick
for k = 1:num_frames
    frame = read(behavior_vid,k);
    frame = frame(:,:,1);
    %display black or white rectangular box on the upper right corner of the
    % image when there is lick or no lick, respectively
    frame(1:60,1:60) = 255*(lick_data(k)==1);
    imshow(frame);
    drawnow;
end



% ANALYZE BEHAVIOR VIDEO
clear all; close all;

% Read maze data
%------------------------------------------------------------
source = 'hand';
lick_source = 'hand-lick.txt';
fps = 20;
%Read the lick file
fileID = fopen(lick_source,'r');
lick_data = fscanf(fileID,'%d');
num_frames = length(lick_data);


% Read video file (kind of slow)
%------------------------------------------------------------
behavior_vid = VideoReader(strcat(source,'.mp4'));
figure,
%the shape to display to indicate lick/no lick
for k = 1:num_frames
    image = read(behavior_vid,k);
    image = image(:,:,1);
    %display black or white rectangular box on the upper right corner of the
    % image when there is lick or no lick, respectively
    image(1:60,1:60) = 255*(lick_data(k)==1);
    imshow(image(:,:,1));
    pause(0.000001);
end



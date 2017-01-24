function make_summary(sources)
%function make_summary(sources)
%
%Combine everything that belongs to a multiday data in a single cell array
%
%   sources : a struct that has the filenames for the maze data, behavior
%   videos, miniscope videos, PCA data, and IC directory. (refer to a 
%   'data_sources.m' file for the specifics of this struct. )
%
%   The function creates a .mat file with the name 'summary+(current 
%   directory name)+.mat'. It contains a cell of two cell arrays:
%     'summary_cell'(of size (#of ICs) x 3) : Columns are for neuron 
%            classification output, ica filter, and center x-y coordinates
%            for each  neuron, respectively.
%     'summary_trial' (of size (#of trials) x 5 ) :  Columns are for start arm,
%            end arm, correctness, ica traces, and spatial trace (extracted 
%            from the behavior video) for each trial, respectively.
%
%Hakan Inan (Jan 15)
%
ic_dir = sources.ic_dir;
ica_source   = get_most_recent_file(ic_dir, 'ica_*.mat');
class_source = get_most_recent_file(ic_dir, 'class_*.txt');
load(ica_source);
class = load_classification(class_source);
num_ics = ica_info.num_ICs;
%%%
%handle the struct 'sources' and extract the number of days from it
%%%
str_sources = fieldnames(sources);
num_fields = length(str_sources);
%find number of occurrences of 'maze*' and 'behavior*'
num_occ_maze = 0;
num_occ_behavior = 0;
for k = 1:num_fields
    if ~isempty(strfind(str_sources{k},'maze'))
        num_occ_maze = num_occ_maze+1;
    end
    if ~isempty(strfind(str_sources{k},'behavior'))
        num_occ_behavior = num_occ_behavior+1;
    end
end
if num_occ_behavior ~= num_occ_maze
    error('Number of maze files needs to match the number of behavior videos');
else
    num_days = num_occ_maze;
end

num_frames = length(imfinfo(sources.miniscope));
size_behavior_vid = [480,640,num_frames];
%%%
%Get trial info from maze output
%%%
%Combine all the days together and treat the whole thing as a single day
trial_start_arms = [];
trial_finish_arms = [];
trial_correctness = [];
compressed_indices = [];
compressed_behavior_vid = zeros(size_behavior_vid,'uint8');
%keep track of the last frame index of the last session to add to the next
last_frame = 0;
%The most time consuming task in the following for loop is compressing the
%behavior video. Notifying the user when it starts and finishes would be
%nice
disp('compressing the behavior video..');
for idx_day = 1:num_days,
    relevant_maze = eval(['sources.maze_',num2str(idx_day)]);
    relevant_behavior_vid_str = eval(['sources.behavior_',num2str(idx_day)]);
    relevant_behavior_vid = VideoReader(relevant_behavior_vid_str); 
    trial_start_arms  = [trial_start_arms;get_trial_start_arm(relevant_maze)];
    trial_finish_arms = [trial_finish_arms;get_trial_finish_arm(relevant_maze)];
    trial_correctness = [trial_correctness;get_trial_correctness(relevant_maze)];
    %compute the trial frames:
    trial_frame_indices = get_trial_frame_indices(relevant_maze, 100, 100);
    compressed_indices = [compressed_indices; last_frame+ ...
        compress_frame_indices(trial_frame_indices)];       
    %read the whole behavior movie:
    whole_behavior_vid = read(relevant_behavior_vid);
    %calculate frames to include in the compressed behavior movie
    num_trials_total = size(trial_frame_indices,1);
    frames = [];
    for k = 1:num_trials_total
        frames = [frames, trial_frame_indices(k,1):trial_frame_indices(k,2)];
    end
    compressed_behavior_vid(:,:,last_frame+1:last_frame+length(frames)) = ...
    squeeze(whole_behavior_vid(:,:,1,frames)); %keep one of the 3 channels
    clear whole_behavior_vid;   
    %update last_frame
    last_frame = last_frame + compressed_indices(end,end);
end
disp('finished compressing the behavior video');
num_trials_total = size(compressed_indices,1);
%%%
%construct 'summary_cell', (#of ICs) x 3 cell array
%%%
%for below 1:cell type 2:cell ica filters 3:cell coordinate
summary_cell = cell(num_ics,3);
sizeImg = size(ica_filters(:,:,1));
%matrices pix_x,pix_y come into play when calculating the center of an IC:
pix_x = repmat(1:sizeImg(2),sizeImg(1),1);
pix_y = repmat(1:sizeImg(1),sizeImg(2),1)';
Aa = zeros(num_ics,2);
for ic_idx = 1:num_ics
    summary_cell{ic_idx,1} = class{ic_idx};
    summary_cell{ic_idx,2} = ica_filters(:,:,ic_idx);
    if (strcmp(class{ic_idx}, 'not a cell'))
        summary_cell{ic_idx,3}=0;summary_cell{ic_idx,3}=0;
        continue;
    end   
    ic_img = ica_filters(:,:,ic_idx);
    ic_img = threshold_ic_filter(ic_img, 0.3);
    %calculate the center of mass of the IC   
    ic_img = ic_img / sum(sum(ic_img)); %make the image sum to 1
    center_x = round(sum(sum(pix_x.*ic_img)));
    center_y = round(sum(sum(pix_y.*ic_img)));
    summary_cell{ic_idx,3} =  [center_x,center_y];
end
%%%
%construct 'summary_trial', (#of trials) x 5 cell array
%%%
ica_traces_t = ica_traces';
%background image to use in the extraction of mouse coordinates
bg_image = imread('bg_image_east.jpg'); 
bg_image = squeeze(bg_image(:,:,1)); % use 1st channel only
%for below, 1:start arm 2:end arm 3:correctness 4:traces 5:spatial trace 
summary_trial = cell(num_trials_total,5);
for trial = 1:num_trials_total
    frames = compressed_indices(trial,:);
    summary_trial{trial,1} = trial_start_arms{trial};
    summary_trial{trial,2} = trial_finish_arms{trial};
    summary_trial{trial,3} = trial_correctness(trial);
    summary_trial{trial,4} = ica_traces_t(:,frames(1):frames(2));
    video = compressed_behavior_vid(:,:,frames(1):frames(2));
    disp(sprintf('getting xy trace for the trial %d',trial));
    summary_trial{trial,5} = get_spatial_trace_from_behavior(bg_image,video);    
end


summary = cell(1,2);
summary{1} = summary_cell;
summary{2} = summary_trial;
%get current folder name:
str = pwd;
for k = length(str):-1:1
    %replace '-' with '_':
    if strcmp(str(k),'-')
        str(k) = '_';
    end
    if strcmp(str(k),'/')||strcmp(str(k),'\')
        folder = str(k+1:end);
        break
    end
end
varname = ['summary_',folder];
eval([varname '=summary;']);
%save summary:
save([varname,'.mat'],varname);



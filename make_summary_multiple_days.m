
clear all; close all;
sources = data_sources;
%Select classified ICs
ic_dir = 'ica001';
ica_source   = get_most_recent_file(ic_dir, 'ica_*.mat');
class_source = get_most_recent_file(ic_dir, 'class_*.txt');
load(ica_source);
class = load_classification(class_source);
num_ics = ica_info.num_ICs;
num_days = str2double(sources.days);
num_frames = length(imfinfo(sources.miniscope));
size_behavior_vid = [480,640,num_frames];
%Get trial info from maze output
%Combine all the days together and treat the whole thing as a single day
sources = data_sources;
trial_start_arms = [];
trial_finish_arms = [];
trial_correctness = [];
compressed_indices = [];
compressed_behavior_vid = zeros(size_behavior_vid,'uint8');
frame_to_start = 0;
for idx_day = 1:num_days,
    relevant_maze = eval(['sources.maze_',num2str(idx_day)]);
    relevant_behavior_vid = eval(['sources.behavior_',num2str(idx_day)]);
    behavior_vid = VideoReader(relevant_behavior_vid); 
    trial_start_arms  = [trial_start_arms;get_trial_start_arm(relevant_maze)];
    trial_finish_arms = [trial_finish_arms;get_trial_finish_arm(relevant_maze)];
    trial_correctness = [trial_correctness;get_trial_correctness(relevant_maze)];
    %compute the trial frames:
    trial_frame_indices = get_trial_frame_indices(relevant_maze, 100, 100);
    compressed_indices = [compressed_indices; frame_to_start+ ...
        compress_frame_indices(trial_frame_indices)];       
    %read the whole behavior movie:
    whole_behavior_vid = read(behavior_vid);
    %calculate frames to include in the compressed behavior movie
    num_trials = size(trial_frame_indices,1);
    frames = [];
    for k = 1:num_trials
        frames = [frames, trial_frame_indices(k,1):trial_frame_indices(k,2)];
    end
    compressed_behavior_vid(:,:,frame_to_start+1:frame_to_start+length(frames)) = ...
    squeeze(whole_behavior_vid(:,:,1,frames));
    clear whole_behavior_vid;   
    frame_to_start = frame_to_start + compressed_indices(end,end);
end


num_trials = size(compressed_indices,1);

%%%
%'summary_cell' is the summary of the whole recording session
%%%
%for below 1:cell type 2:cell ica filters 3:cell coordinate
summary_cell = cell(num_ics,3);
sizeImg = size(ica_filters(:,:,1));
%matrices gravx,gravy come into play when calculating the center of an IC:
gravx = repmat(1:sizeImg(2),sizeImg(1),1);
gravy = repmat(1:sizeImg(1),sizeImg(2),1)';
Aa = zeros(num_ics,2);
for ic_idx = 1:num_ics
    summary_cell{ic_idx,1} = class{ic_idx};
    summary_cell{ic_idx,2} = ica_filters(:,:,ic_idx);
    if (strcmp(class{ic_idx}, 'not a cell'))
        summary_cell{ic_idx,3}=0;summary_cell{ic_idx,3}=0;
        continue;
    end   
    %%%
    %extracting coordinates:
    %%%
    a = ica_filters(:,:,ic_idx);
    a = threshold_ic_filter(a, 0.3);
    %calculate the center of gravity of the IC in the image:  
    normalizer = sum(sum(a));
    centerx = round(sum(sum(gravx.*a))/normalizer);
    centery = round(sum(sum(gravy.*a))/normalizer);
    summary_cell{ic_idx,3} =  [centerx,centery];
end
%%%
%traces and trial labels for each cell:
%%%
wholeSession = ica_traces';
east_img = imread('bg_image_east.jpg');
east_img = squeeze(east_img(:,:,1));
west_img = imread('bg_image_west.jpg');
west_img = squeeze(west_img(:,:,1));
%for below, 1:start arm 2:end arm 3:correctness 4:traces 5:spatial trace 
summary_trial = cell(num_trials,5);
for trial = 1:num_trials
    frames = compressed_indices(trial,:);
    summary_trial{trial,1} = trial_start_arms{trial};
    summary_trial{trial,2} = trial_finish_arms{trial};
    summary_trial{trial,3} = trial_correctness(trial);
    summary_trial{trial,4} = wholeSession(:,frames(1):frames(2));
    video = compressed_behavior_vid(:,:,frames(1):frames(2));
    if strcmp('east',trial_start_arms{trial})
        bg_image = east_img;
    else
        bg_image = east_img;
    end 
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



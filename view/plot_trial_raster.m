function [cell_order] = plot_trial_raster(traces,frames,trial_num,idx)

% plot_trial_raster plots all normalized cell traces over the course of a
% single trial
%
% Inputs:
%   traces: [num_frames (over all trials) x num_cells] 
%       (matrix of traces for each cell over all concatenated trials)
%   frames: [num_trials x 4]
%   trial_num = index of trial to plot
%   varagin = [1 x numCells]
%       (matrix containing the cell order to use in the raster plot)

for i=1:size(traces,2)
    trial_frames = traces(frames(trial_num,1):frames(trial_num,end),i);
    frames_norm = trial_frames./abs(max(trial_frames));
    norm_trial(i,:) = frames_norm.';
    clear frames_norm trial_frames
end

if(~exist('idx','var'))
    [M,I] = max(norm_trial,[],2);
    [B,idx] = sort(I);
end

cell_order = idx;

for i=1:length(idx)
    new_norm_trial(i,:) = norm_trial(idx(i),:);
end

new_norm_trial(:,frames(trial_num,2)-frames(trial_num,1)) = 1;
new_norm_trial(:,frames(trial_num,3)-frames(trial_num,1)) = 1;

figure
imagesc(new_norm_trial,[0 1]);
colorbar
xlabel('Frame')
ylabel('Cell');
str = sprintf('Trial %i',trial_num);
title(str);
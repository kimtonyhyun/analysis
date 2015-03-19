function [cellOrder] = plot_trial_raster(ic_traces,binned_frames,trial_num)

for i=1:size(ic_traces,2)
    frames = ic_traces(binned_frames(trial_num,1):binned_frames(trial_num,end),i);
    frames_norm = frames./abs(max(frames));
    norm_trial(i,:) = frames_norm.';
    clear frames_norm frames
end

[M,I] = max(norm_trial,[],2);
[B,idxN] = sort(I);
cellOrder = idxN;

for i=1:length(I)
    new_norm_trial(i,:) = norm_trial(idxN(i),:);
end

new_norm_trial(:,binned_frames(trial_num,2)-binned_frames(trial_num,1)) = 1;
new_norm_trial(:,binned_frames(trial_num,3)-binned_frames(trial_num,1)) = 1;

figure(2)
imagesc(new_norm_trial,[0 1]);
colorbar
xlabel('Frame')
ylabel('Cell');
str = sprintf('Trial %i',trial_num);
title(str);
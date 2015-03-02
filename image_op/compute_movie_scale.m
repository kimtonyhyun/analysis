function clim = compute_movie_scale(M)
% Compute an appropriate viewing scale (CLim) for the provided movie

[height, width, ~] = size(M);
maxVec = single(reshape(max(M,[],3), height*width, 1));
minVec = single(reshape(min(M,[],3), height*width, 1));
quantsMax = quantile(maxVec,[0.85,0.87,0.9,0.93,0.95]);
quantsMin = quantile(minVec,[0.85,0.87,0.9,0.93,0.95]);
clim = [mean(quantsMin),mean(quantsMax)]*1.1;
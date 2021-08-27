function clim = compute_movie_scale(M)
% Compute an appropriate viewing scale (CLim) for the provided movie

switch ndims(M)
    case 3 % Standard movie
        [height, width, ~] = size(M);
        
        maxVec = reshape(max(M,[],3), height*width, 1);
        minVec = reshape(min(M,[],3), height*width, 1);
        
        quantsMax = quantile(maxVec,[0.85,0.87,0.9,0.93,0.95]);
        quantsMin = quantile(minVec,[0.85,0.87,0.9,0.93,0.95]);

    case 2 % Image
        quantsMax = quantile(M(:),[0.95, 0.99]);
        quantsMin = quantile(M(:),[0.01, 0.05]);
end

clim = [mean(quantsMin) mean(quantsMax)];
clim_range = clim(2)-clim(1);
clim = clim + 0.1*clim_range*[-1 1];

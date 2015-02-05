function M_norm = norm_by_background_mean(movie)
%Normalizes every frame in the movie by the mean pixel value of the
%background to mitigate the temporal variation in the brightness
%
%   movie: movie matrix , [h x w x num_frames]
%
% Hakan Inan, 15-Jan-4
%

M_norm = zeros(size(movie),'single');
[height,width,numFrames] = size(movie);

for k = 1:numFrames
    if (mod(k,1000)==0)
        fprintf('  Frames %d of %d done\n', k, numFrames);
    end
    Frame = reshape(movie(:,:,k),height*width,1);
    thresh = quantile(Frame,0.7);
    Z = mean(Frame(Frame<thresh));
    M_norm(:,:,k) = movie(:,:,k)/Z; %normalize by Z
end


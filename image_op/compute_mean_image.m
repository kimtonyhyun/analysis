function A = compute_mean_image(M)
% Compute the mean projection image of the movie 'M'. The movie 'M' can be
% the name of a file.

if ~ischar(M)
    [height, width, num_frames] = size(M);
    A = zeros(height, width);
    
    for k = 1:num_frames
        A = A + M(:,:,k);
    end
    A = A / num_frames;
else % HDF5 filename (string)
    [~, A] = compute_std_image(M);
end


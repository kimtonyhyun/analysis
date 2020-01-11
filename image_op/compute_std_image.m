function [S, A] = compute_std_image(M)
% Compute the per-pixel standard deviation in a way that is more memory
% efficient than S = std(M,[],3);
%
% The variable 'M' can be a movie in memory or a path to a HDF5 file

if ~ischar(M)
    [height, width, num_frames] = size(M);
    A = zeros(height, width);
    A2 = zeros(height, width);

    for k = 1:num_frames
        A = A + M(:,:,k);
        A2 = A2 + M(:,:,k).^2;
    end
else % HDF5 filename (string)
    % Default dataset name for the movie
    movie_dataset = '/Data/Images';

    % Grab the movie parameters
    [movie_size, ~] = get_dataset_info(M, movie_dataset);
    height = movie_size(1);
    width = movie_size(2);
    num_frames = movie_size(3);
    
    A = zeros(height, width);
    A2 = zeros(height, width);
    
    frame_chunk_size = 2500;
    [frame_chunks, num_chunks] = make_frame_chunks(num_frames, frame_chunk_size);
    
    for i = 1:num_chunks
        fprintf('%s: Reading frames %d to %d (out of %d)...\n',...
            datestr(now), frame_chunks(i,1), frame_chunks(i,2), num_frames);

        chunk_start = frame_chunks(i,1);
        chunk_count = frame_chunks(i,2) - frame_chunks(i,1) + 1;

        movie_chunk = h5read(M, movie_dataset,...
                             [1 1 chunk_start],...
                             [height width chunk_count]);

        for k = 1:size(movie_chunk,3)
            A = A + movie_chunk(:,:,k);
            A2 = A2 + movie_chunk(:,:,k).^2;
        end
    end
end

A = A / num_frames;
A2 = A2 / num_frames;
S = sqrt(A2-A.^2);
function save_movie_to_hdf5(movie, outname, dataset_name)
% Save a movie matrix [height x width x num_frames] into a HDF5 file.
%
% 2015 01 28 Tony Hyun Kim

[height, width, num_frames] = size(movie);
data_type = class(movie);

chunk_size = [height width 1];

h5create(outname, dataset_name, [height width num_frames],...
    'DataType', data_type,...
    'ChunkSize', chunk_size);
h5write(outname, dataset_name, movie);

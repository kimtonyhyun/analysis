function M = load_movie_from_hdf5(source, dataset_name)
% Load a movie matrix [height x width x num_frames] from a HDF5 file.
%
% 2015 01 28 Tony Hyun Kim

M = h5read(source, dataset_name)

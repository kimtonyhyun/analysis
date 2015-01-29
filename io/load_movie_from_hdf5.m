function M = load_movie_from_hdf5(source, varargin)
% Load a movie matrix [height x width x num_frames] from a HDF5 file.
%   Optional parameter to specify the dataset name.
%
% 2015 01 28 Tony Hyun Kim
if isempty(varargin)
    dataset_name = '/Data/Images';
else
    dataset_name = varargin{1};
end

M = h5read(source, dataset_name);

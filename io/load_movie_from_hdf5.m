function M = load_movie_from_hdf5(source, varargin)
% Load a movie matrix [height x width x num_frames] from a HDF5 file.
%   Optional parameter to specify the frame range as [start end]
%

movie_dataset = '/Data/Images';

if ~isempty(varargin) % Read a subset of the movie
    frame_range = varargin{1}; % [Start End]
    frame_start = frame_range(1);
    frame_count = frame_range(2) - frame_range(1) + 1;
    
    [movie_size, ~] = get_dataset_info(source, movie_dataset);
    
    M = h5read(source, movie_dataset,...
               [1 1 frame_start],...
               [movie_size(1) movie_size(2) frame_count]);
    
else % Read the whole file
    M = h5read(source, movie_dataset);
end
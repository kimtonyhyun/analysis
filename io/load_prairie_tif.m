function M = load_prairie_tif(path_to_tif, varargin)
% Prairie saves frames of a series as separate TIF files, where each file
% contains exactly one image. This loader simply enumerates all TIF files
% in a given directory, and return the result as a single movie matrix.
% Optional arguments can be used to specify loading specific channels or
% specific planes of a multi-plane dataset.
%
% Usage:
%   M = load_prairie_tif(pwd);
%
% TODO:
%   Consider parsing XML file, if available.

channel_ind = [];
slice_ind = [];
for k = 1:length(varargin)
    vararg = varargin{k};
    if ischar(vararg)
        switch lower(vararg)
            case {'ch', 'channel'}
                channel_ind = varargin{k+1};
            case {'sl', 'slice'} % Used in volumetric data
                slice_ind = varargin{k+1};
        end
    end
end

% Build file pattern
file_pattern = '*.tif'; % Default: All TIF files in directory
if (~isempty(channel_ind) && ~isempty(slice_ind))
    file_pattern = sprintf('*_Ch%d_%06d.*.tif', channel_ind, slice_ind);
elseif (~isempty(channel_ind))
    file_pattern = sprintf('*_Ch%d_*.tif', channel_ind);
elseif (~isempty(slice_ind))
    file_pattern = sprintf('*_%06d.*.tif', slice_ind);
end
fprintf('%s: Looking for files matching "%s" in directory "%s"... ',...
    datestr(now), file_pattern, path_to_tif);

files = dir(fullfile(path_to_tif, file_pattern));
num_frames = length(files);

fprintf('Found %d images!\n', num_frames);

% Open one image to get the X, Y dimensions and pixel data type
A = imread(fullfile(path_to_tif, files(1).name));
M = zeros(size(A,1), size(A,2), num_frames, class(A));

for i = 1:num_frames
    M(:,:,i) = imread(fullfile(path_to_tif, files(i).name));
end
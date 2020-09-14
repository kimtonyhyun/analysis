function M = load_movie(source)
% Call the appropriate movie loading method based on the file extension
%   of the specified `source`
%
% 2015 01 28 Tony Hyun Kim

[~, ~, ext] = fileparts(source);
ext = ext(2:end); % Remove the leading dot

if isempty(ext)
    files = dir(strcat(source, '.*'));
    if isempty(files)
        error('File not found!');
    else
        % TODO: What to do if more than 1 match?
        source = files(1).name;
        [~, ~, ext] = fileparts(source);
        ext = ext(2:end);
    end
end

switch lower(ext)
    case {'tif', 'tiff'}
        M = load_movie_from_tif(source);
    case {'h5', 'hdf5'}
        M = load_movie_from_hdf5(source);
    case {'dcimg'}
        M = load_movie_from_dcimg(source);
    otherwise
        error('Unrecognized file type!');
end
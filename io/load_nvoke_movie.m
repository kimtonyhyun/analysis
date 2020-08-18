function movie = load_nvoke_movie(source)
% Based on 'load_movie_from_tif', which was originally written to work with
% nVista v2 TIF outputs (~2015-2016).
%
% Adapted to load nVoke v1 HDF5/XML outputs (~2017). The use case (c.f.
% direct 'h5read') is to fill in dropped frames.

% nVoke HDF5 stores image data in '/images'
movie = h5read(source, '/images');
[height, width, num_hdf5_frames] = size(movie);

% Parse XML
xml_filename = convert_extension(source, 'xml');
xml_struct = parse_miniscope_xml(xml_filename);

% Miniscope XML file tallies recorded and dropped frames separately.
% Make sure that the recorded frames match what is in the HDF5 file.
num_recorded_frames = str2double(xml_struct.frames);
assert(num_recorded_frames == num_hdf5_frames,...
    '  Unexpected number of frames in HDF5 file!');

num_dropped_frames = str2double(xml_struct.dropped_count);
fprintf('Found %d dropped frames in "%s"\n', num_dropped_frames, source);
if (num_dropped_frames > 0)
    dropped_frames = str2num(xml_struct.dropped); %#ok<ST2NM>

    num_total_frames = num_recorded_frames + num_dropped_frames;
    movie_corr = zeros(height, width, num_total_frames, 'like', movie);

    % Each missing frame slot will be replaced with the PREVIOUS
    % recorded frame. Except when the first frames of a recording are
    % dropped; in that case, we replace with the SUBSEQUENT recorded
    % frame.
    idx = 0;
    for k = 1:num_total_frames
        if ismember(k, dropped_frames) % Is a dropped frame
            sample_idx = max(1, idx); % Special handling for first frame
            movie_corr(:,:,k) = movie(:,:,sample_idx);
        else
            idx = idx + 1;
            movie_corr(:,:,k) = movie(:,:,idx);
        end
    end
    assert(idx == num_recorded_frames,...
        '  Not all recorded frames have been transferred to dropped-frame corrected movie!');

    movie = movie_corr;
end
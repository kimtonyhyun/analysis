function crop_movie(movie_in, movie_out)
% Crop HDF5 movie from file ('movie_in') to file ('movie_out'). The
% cropping region is defined by the user interactively.
%
% If 'movie_out' is left as an empty string, then default name will be
% provided.
%
% The cropping parameters will also be saved in the output HDF5 file under
% the '/Crop' directory
%
% Inputs:
%   movie_in:  Name of incoming HDF5 movie
%   movie_out: Name of outgoing HDF5 movie
%
% Example usage:
%   crop_movie('c9m7d12.hdf5','');
%

if isempty(movie_out)
    [~, name] = fileparts(movie_in);
    movie_out = sprintf('%s_cr.hdf5', name);
    fprintf('crop_movie: Output movie will be saved as "%s"\n', movie_out);
end

% Default dataset name for the movie
movie_dataset = '/Data/Images';

% Grab the movie parameters
h5att = h5info(movie_in, movie_dataset);
in_data_type = datatype_hdf5_to_matlab(h5att.Datatype.Type); % e.g. 'uint16'

height = h5att.Dataspace.Size(1);
width = h5att.Dataspace.Size(2);
num_frames = h5att.Dataspace.Size(3);

% Let the user select the crop region
% TODO: Allow the cropping window to be defined as a parameter to the function
%------------------------------------------------------------
ref_idx = 1;
ref_frame = h5read(movie_in, movie_dataset, [1 1 ref_idx], [height width 1]);
imagesc(ref_frame);
colormap gray;
axis image;
title(strrep(movie_in, '_', '\_'));
fprintf('crop_movie: Please provide a rectangular region over the image.\n');
fprintf('  Double click on the rectangle when done.\n');
h_rect = imrect;
rect_params = round(wait(h_rect));

% Show the cropped image
x0 = rect_params(1); x1 = rect_params(1)+rect_params(3)-1;
y0 = rect_params(2); y1 = rect_params(2)+rect_params(4)-1;
ref_frame_cropped = ref_frame(y0:y1, x0:x1);
imagesc(ref_frame_cropped);
axis image;
title(sprintf('%s (cropped)',strrep(movie_in, '_', '\_')));

input('crop_movie: Press enter to apply crop to entire movie >> ');

% Prepare output movie
%------------------------------------------------------------
cropped_height = rect_params(4);
cropped_width = rect_params(3);
h5create(movie_out, movie_dataset,...
         [cropped_height cropped_width, num_frames],...
         'ChunkSize', [cropped_height, cropped_width 1],...
         'Datatype', in_data_type); % Preserve data type

copy_hdf5_params(movie_in, movie_out);

h5create(movie_out, '/Crop/CropRect', [1 4], 'Datatype', 'uint16');
h5write(movie_out, '/Crop/CropRect', uint16(rect_params));

% Process the remainder of the movie
frame_chunk_size = 2500;
[frame_chunks, num_chunks] = make_frame_chunks(num_frames, frame_chunk_size);

for i = 1:num_chunks
    fprintf('%s: Cropping frames %d to %d (out of %d)...\n',...
        datestr(now), frame_chunks(i,1), frame_chunks(i,2), num_frames);
    
    chunk_start = frame_chunks(i,1);
    chunk_count = frame_chunks(i,2) - frame_chunks(i,1) + 1;
    
    movie_chunk = h5read(movie_in, movie_dataset,...
                         [1 1 chunk_start],...
                         [height width chunk_count]);
                         
    h5write(movie_out, movie_dataset,...
            movie_chunk(y0:y1, x0:x1, :),...
            [1 1 chunk_start],...
            [cropped_height cropped_width chunk_count]);
end
fprintf('%s: Done!\n', datestr(now));
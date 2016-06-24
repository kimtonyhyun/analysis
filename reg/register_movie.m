function register_movie(movie_in, movie_out, varargin)
% Motion correct HDF5 movie from file ('movie_in') to file ('movie_out').
% Registration is performed with TurboReg (turbocoreg.mex). The
% registration is performed after filtering the frames with a filter.
% 
% The default filter (called "mosaic") has worked well for Miniscope and
% other 1p movies. The Mosaic filter is the default option.
%
% If 'movie_out' is left as an empty string, then default name will be
% provided.
%
% Inputs:
%   movie_in:  Name of incoming HDF5 movie
%   movie_out: Name of outgoing HDF5 movie
%
% Example usage:
%   register_movie('c9m7d12.hdf5', '');
%

if isempty(movie_out)
    [~, name] = fileparts(movie_in);
    movie_out = sprintf('%s_mc.hdf5', name);
    fprintf('register_movie: Output movie will be saved as "%s"\n', movie_out);
end

% Default dataset name for the movie
movie_dataset = '/Data/Images';

% Grab the movie parameters
[movie_size, ~] = get_dataset_info(movie_in, movie_dataset);
height = movie_size(1);
width = movie_size(2);
num_frames = movie_size(3);

% Begin TurboReg processing
%------------------------------------------------------------
ref_idx = 1; % By default, movie is registered against the first frame
im_ref = h5read(movie_in, movie_dataset, [1 1 ref_idx], [height width 1]);

filter_type = 'mosaic'; % Default filter option

for k = 1:length(varargin)
    vararg = varargin{k};
    if ischar(vararg)
        vararg = lower(vararg);
        switch vararg
            case 'ref' % Externally provided reference image
                im_ref = varargin{k+1};
            otherwise
                filter_type = vararg;
        end
    end
end

switch filter_type
    case 'mosaic'
        ssm_radius = 20;
        asm_radius = 5;
        hDisk  = fspecial('disk', ssm_radius);
        hDisk2 = fspecial('disk', asm_radius);
        transform = @(A) mosaic_transform(A, hDisk, hDisk2);

    case {'identity', 'nofilter'}
        transform = @(A) identity_transform(A);
        
    otherwise
        fprintf('  Error: Did not recognize filter type "%s"!\n', filter_type);
        return
end
fprintf('register_movie: Using "%s" filter...\n', filter_type);
im_ref = transform(single(im_ref));

% Specify ROI
fprintf('register_movie: Please select ROI for TurboReg\n');
imagesc(im_ref); axis image; colormap gray;
title(strrep(movie_in,'_','\_'));
h_poly = impoly;
mask_xy = getPosition(h_poly);
mask_ref = single(poly2mask(mask_xy(:,1), mask_xy(:,2), height, width));

% Turboreg options
options.rotation_enable = false;
options.mingain = 0.0; % Max accuracy
options.levels = calculate_pyramid_depth(min(height, width));

input('register_movie: Please enter to proceed >> ');

% Prepare output movie
%------------------------------------------------------------
h5create(movie_out, movie_dataset,...
         [height width num_frames],...
         'ChunkSize', [height width 1],...
         'Datatype', 'single');
     
copy_hdf5_params(movie_in, movie_out);     
     
h5create(movie_out, '/MotCorr/MaskXY', size(mask_xy), 'Datatype', 'double');
h5write(movie_out, '/MotCorr/MaskXY', mask_xy);

% Apply TurboReg
frame_chunk_size = 500;
[frame_chunks, num_chunks] = make_frame_chunks(num_frames, frame_chunk_size);

for i = 1:num_chunks
    fprintf('%s: Registering frames %d to %d (out of %d)...\n',...
        datestr(now), frame_chunks(i,1), frame_chunks(i,2), num_frames);
    
    chunk_start = frame_chunks(i,1);
    chunk_count = frame_chunks(i,2) - frame_chunks(i,1) + 1;
    
    movie_chunk = h5read(movie_in, movie_dataset,...
                         [1 1 chunk_start],...
                         [height width chunk_count]);

	movie_chunk_reg = zeros(size(movie_chunk), 'single');
    
    for frame_idx = 1:size(movie_chunk,3);
        im_coreg = single(movie_chunk(:,:,frame_idx));
        im_reg = transform(im_coreg);
        
        movie_chunk_reg(:,:,frame_idx) = ...
            turbocoreg(im_ref, mask_ref, im_reg, im_coreg, options);
    end
    
    h5write(movie_out, movie_dataset,...
            movie_chunk_reg,...
            [1 1 chunk_start],...
            [height width chunk_count]);
end
fprintf('%s: Done!\n', datestr(now));

end % register_movie

function depth = calculate_pyramid_depth(len)
    min_size = 45;
    depth = 0;
    while (min_size <= len)
        len = len/2;
        depth = depth + 1;
    end
end

%------------------------------------------------------------
% Built-in image filters
%------------------------------------------------------------
function A_tr = mosaic_transform(A, ssm_filter, asm_filter)
    % Derived from Inscopix's Mosaic. Frame is transformed in two steps:
    %   (1) "Subtract spatial mean", then
    %   (2) "Apply spatial mean"
    A_tr = A - imfilter(A, ssm_filter, 'replicate');
    A_tr = imfilter(A_tr, asm_filter);
end

function A_tr = identity_transform(A)
    A_tr = A;
end
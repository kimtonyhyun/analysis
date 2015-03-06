function view_movie(M, varargin)
% Displays the frames of a movie matrix M [height x row x num_frames]
%   (Note: also works with a single image)

num_repeats = 1;
rescale_each_frame = false;
if ~isempty(varargin)
    len = length(varargin);
    for k = 1:len
        switch lower(varargin{k})
            case 'repeats'
                num_repeats = varargin{k+1};
            case 'rescale'
                rescale_each_frame = true;
        end
    end
end

num_frames = size(M,3);

if (rescale_each_frame || isa(M, 'uint16'))
    % Every frame is rescaled for the raw movie
    h = imagesc(M(:,:,1));
else % Otherwise, use common CLim scaling
    movie_clim = compute_movie_scale(M);
    h = imagesc(M(:,:,1), movie_clim);
end
axis image;
truesize;
colormap gray;
xlabel('x [px]');
ylabel('y [px]');

for r = 1:num_repeats
    for k = 1:num_frames
        title(sprintf('Frame %d of %d', k, num_frames));
        set(h, 'CData', M(:,:,k));
        drawnow;
    end
end

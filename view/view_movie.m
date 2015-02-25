function view_movie(M, varargin)
% Displays the frames of a movie matrix M [height x row x num_frames]
%   (Note: also works with a single image)
% Optional input will repeat the movie.
if isempty(varargin)
    num_repeats = 1;
else
    num_repeats = varargin{1};
end

num_frames = size(M,3);

movie_clim = compute_movie_scale(M);
h = imagesc(M(:,:,1), movie_clim);
axis image;
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

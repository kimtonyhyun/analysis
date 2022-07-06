function view_movie(M, varargin)
% Displays the frames of a movie matrix M. Movie can be:
%   - Grayscale: [height x width x num_frames]
%   - RGB: [height x width x RGB x num_frames]

movie_clim = [];
use_outline = 0;
poi = [];
num_repeats = inf;

for k = 1:length(varargin)
    vararg = varargin{k};
    if ischar(vararg)
        switch lower(vararg)
            case {'repeat', 'repeats'}
                num_repeats = varargin{k+1};
            case 'poi'
                poi = varargin{k+1};
            case 'clim'
                movie_clim = varargin{k+1};
            case 'boundary'
                use_outline = 1;
                ds = varargin{k+1};
        end
    end
end

if isempty(movie_clim)
    switch class(M)
        case {'uint8'}
            movie_clim = [0 255];
        case {'uint16'}
            movie_clim = [0 0.9*max(M(:))];
        case {'single', 'double'}
            movie_clim = compute_movie_scale(M);
    end
end

nd = ndims(M);
figure;
switch nd
    case 3 % [height x width x num_frames]       
        get_frame = @(k) M(:,:,k);
        h = imagesc(get_frame(1), movie_clim);
        colormap gray;
    case 4 % [height x width x RGB x num_frames]
        get_frame = @(k) M(:,:,:,k);
        h = image(get_frame(1));
end
ht = get(gca, 'Title');
num_frames = size(M,nd);

axis image;
truesize;
xlabel('x [px]');
ylabel('y [px]');

% Toolbar implementation in Matlab 2018b+ is broken
if ~verLessThan('matlab', '9.5')
    addToolbarExplorationButtons(gcf);
    ax = gca;
    set(ax.Toolbar, 'Visible', 'off');
end

if use_outline
    hold on;
    plot_boundaries(ds, 'Color', 'g');
    hold off;
end

if ~isempty(poi)
    hold on;
    plot(poi(:,1), poi(:,2), 'm*');
    hold off;
end

num_playbacks = 1;
while (num_playbacks <= num_repeats) 
    for k = 1:num_frames
        try
            set(h, 'CData', get_frame(k));
            ht.String = sprintf('Frame %d of %d', k, num_frames);
        catch
            cprintf('blue', 'view_movie terminated by user\n');
            return;
        end
        drawnow;
    end
    num_playbacks = num_playbacks + 1;
end

fprintf('test');
function run_pca(movie_source, num_PCs)

fprintf('%s: Loading %s...\n', datestr(now), movie_source);
M = load_movie(movie_source);
[height, width, num_frames] = size(M);

% Make each frame zero-mean in place
fprintf('%s: Frame-by-frame normalization of movie...\n', datestr(now));
F = compute_mean_fluorescence(M);
F = reshape(F,1,1,num_frames);
M = bsxfun(@minus, M, F);

% Reshape movie into [space x time] matrix
num_pixels = height * width;
M = reshape(M, num_pixels, num_frames);

% PCA
%------------------------------------------------------------
[pca_filters, pca_traces, S] = compute_pca(M, num_PCs); %#ok<*NASGU,*ASGLU>

savename = sprintf('pca_n%d.mat', num_PCs);

pca_info.movie_height = height;
pca_info.movie_width  = width;
pca_info.movie_frames = num_frames;
pca_info.num_PCs = num_PCs; %#ok<STRNU>
save(savename, 'pca_info', 'pca_filters', 'pca_traces', 'S');
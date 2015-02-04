function run_pca(movie_source, num_PCs)

fprintf('%s: Loading %s...\n', datestr(now), movie_source);
M = load_movie(movie_source);
[height, width, num_frames] = size(M);

% Reshape movie into [space x time] matrix
num_pixels = height * width;
M = reshape(M, num_pixels, num_frames);

% Make each frame zero-mean in place
fprintf('%s: Frame-by-frame normalization of movie...\n', datestr(now));
for i = 1:num_frames
    M(:,i) = M(:,i) - mean(M(:,i));
end

% PCA
%------------------------------------------------------------
[pca_filters, pca_traces, S] = compute_pca(M, num_PCs); %#ok<*NASGU,*ASGLU>

savename = sprintf('pca_n%d.mat', num_PCs);

pca_info.movie_height = height;
pca_info.movie_width  = width;
pca_info.movie_frames = num_frames;
pca_info.num_PCs = num_PCs; %#ok<STRNU>
save(savename, 'pca_info', 'pca_filters', 'pca_traces', 'S');

fprintf('%s: All done!\n', datestr(now));
function run_pca(movie, num_PCs)

[height, width, num_frames] = size(movie);

% PCA
%------------------------------------------------------------
[pca_filters, pca_traces, S] = compute_pca(movie, num_PCs); %#ok<*NASGU,*ASGLU>

savename = sprintf('pca_n%d.mat', num_PCs);

pca_info.movie_height = height;
pca_info.movie_width  = width;
pca_info.movie_frames = num_frames;
pca_info.num_PCs = num_PCs; %#ok<STRNU>
save(savename, 'pca_info', 'pca_filters', 'pca_traces', 'S');
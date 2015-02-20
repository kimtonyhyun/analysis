function [ica_filters, ica_traces] = compute_ica_pairs(pca_source, icaw_source)

% Load PCA results
load(pca_source, 'pca_info', 'pca_filters', 'pca_traces');
height  = pca_info.movie_height;
width   = pca_info.movie_width;

% Load ICA weights
load(icaw_source, 'ica_W');
num_ICs = size(ica_W, 1);

% Compute ICA pairs
ica_traces  = ica_W*pca_traces;
ica_filters = ica_W*pca_filters;

% Format the ICA dimensions for output
ica_traces  = ica_traces'; % [time x IC]
ica_filters = reshape(ica_filters, num_ICs, height, width);
ica_filters = permute(ica_filters, [2 3 1]); % [height, width, IC]

% Sort ICs based on skewness
skews = skewness(ica_traces);
[~, sorted] = sort(skews, 'descend');
ica_traces  = ica_traces(:,sorted); %#ok<*NASGU>
ica_filters = ica_filters(:,:,sorted);
clear skews sorted;
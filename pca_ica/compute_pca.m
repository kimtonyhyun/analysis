function [pca_filters, pca_traces, S] = compute_pca(movie, num_PCs)

% For understanding the logic behind these operations, the Wikipedia
%   article on Singular Value Decomposition is recommended

[height, width, num_frames] = size(movie);
num_pixels = height * width;

% Reshape movie into [space x time] matrix
movie = reshape(movie, num_pixels, num_frames);

% Compute the covariance matrix [time x time]
fprintf('%s: Computing covariance matrix...\n', datestr(now));
C = cov(movie, 1);    % Normalized by num_pixels
C = num_pixels*C; % Undo the normalization

fprintf('%s: Computing temporal PCs...\n', datestr(now));
options.issym = 'true';
C = double(C); % Conversion needed for 'eigs'
[pca_traces, cov_eigs] = eigs(C, num_PCs, 'LM', options);

cov_eigs = diag(cov_eigs)'; % Don't need the matrix

% Keep only positive eigenvalues. Just a safeguard, should not have
% non-positive eigenvalues, unless the covariance matrix is pathological
sieve = cov_eigs > 0;
pca_traces = pca_traces(:, sieve);
cov_eigs = cov_eigs(:, sieve);
num_PCs = sum(sieve);
clear sieve;

% Singular values
S = diag(cov_eigs.^(1/2));

% Compute the corresponding spatial PCs
fprintf('%s: Computing corresponding PC filters...\n', datestr(now));
movie = movie - repmat(mean(movie,1), num_pixels, 1); % Space normalized
pca_filters = (movie * pca_traces) / S;

% Perform explicit normalization of the PCA filters
for pc_idx = 1:num_PCs
    pca_filter = pca_filters(:, pc_idx);
    pca_filters(:,pc_idx) = pca_filter / norm(pca_filter);
end

% Output formatting
pca_traces  = single(pca_traces);
pca_traces  = pca_traces';  % [num_PCs x time]
pca_filters = pca_filters'; % [num_PCs x space]

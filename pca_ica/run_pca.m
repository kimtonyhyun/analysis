function run_pca(movie_source, num_PCs,varargin)
% Runs PCA factorization of the movie provided in `movie_source`. Saves the
% result to 'pca_(...).mat' file.
%
% Inputs:
%   movie_source: Movie filename
%   num_PCs: Number of requested PCs
%
% Variable input arguments:
%   'trim': if added as an argument then the script trims pixels in the 
%       movie that do not cross above the median of the maximum pixel 
%       values.
%
% Example usage:
%   run_pca('c9m7d25_dff.hdf5', 500,'trim');
%

% Defaults
do_trim = 0;

if ~isempty(varargin)
    len = length(varargin);
    for k = 1:len
        switch varargin{k}
            case 'trim'
                do_trim = 1;
        end
    end
end


fprintf('%s: Loading %s...\n', datestr(now), movie_source);
M = load_movie(movie_source);
[height, width, num_frames] = size(M);

% Reshape movie into [space x time] matrix
num_pixels = height * width;
M = reshape(M, num_pixels, num_frames);

% Make each frame zero-mean in place
mean_M = mean(M,1);
M = bsxfun(@minus, M, mean_M);

if do_trim
    max_proj = max(M,[],2);
    idx_kept = find(max_proj>median(max_proj));
    M = M(idx_kept,:);  %#ok<FNDSB>
end

% PCA
%------------------------------------------------------------
[pca_filters, pca_traces, S] = compute_pca(M, num_PCs); %#ok<*NASGU,*ASGLU>

savename = sprintf('pca_n%d.mat', num_PCs);

pca_info.movie_height = height;
pca_info.movie_width  = width;
pca_info.movie_frames = num_frames;
pca_info.num_PCs = num_PCs; 
pca_info.trimmed = do_trim; %#ok<STRNU>

if do_trim
    save(savename,'pca_info', 'pca_filters', 'pca_traces', 'S','idx_kept');
else
    save(savename, 'pca_info', 'pca_filters', 'pca_traces', 'S');
end

fprintf('%s: All done!\n', datestr(now));
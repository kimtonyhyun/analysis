function run_pca(movie_source, num_PCs, varargin)
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
%   'medfilt': Add it as an argument to perform median-filtering on the 
%       movie on a per-frame basis before PCA. 
%
% Example usage:
%   run_pca('c9m7d25_dff.hdf5', 500,'trim','medfilt');
%

% Defaults
do_trim = 0;
do_medfilt = 0;
medfilt_halfwidth = 1;

if ~isempty(varargin)
    for k = 1:length(varargin)
        switch varargin{k}
            case 'trim'
                do_trim = 1;
            case 'medfilt'
                do_medfilt = 1;
        end
    end
end


fprintf('%s: Loading %s...\n', datestr(now), movie_source);
M = load_movie(movie_source);
[height, width, num_frames] = size(M);

% Median filter the movie (optional)
if do_medfilt
    medfilt_neighborhood = (1+2*medfilt_halfwidth)*[1 1];

    for idx_frame = 1:num_frames
        frame = M(:,:,idx_frame);
        M(:,:,idx_frame) = medfilt2(frame, medfilt_neighborhood);
        if mod(idx_frame,1000)== 0
            fprintf('%s: Median-filtered %d frames (out of %d)...\n',...
                datestr(now),idx_frame, num_frames);
        end
    end

    fprintf('%s: Finished median filtering!\n', datestr(now));
end

% Reshape movie into [space x time] matrix
num_pixels = height * width;
M = reshape(M, num_pixels, num_frames);

% Make each frame zero-mean in place
mean_M = mean(M,1);
M = bsxfun(@minus, M, mean_M);

idx_kept = 1:num_pixels;
if do_trim
    max_proj = max(M,[],2);
    idx_kept = find(max_proj>median(max_proj));
    M = M(idx_kept,:);
end

% PCA
%------------------------------------------------------------
[pca_filters, pca_traces, S] = compute_pca(M, num_PCs); %#ok<*NASGU,*ASGLU>
S = diag(S); % Save only the diagonal of S

savename = sprintf('pca_n%d.mat', num_PCs);

pca_info.movie_height = height;
pca_info.movie_width  = width;
pca_info.movie_frames = num_frames;
pca_info.num_PCs = num_PCs; 

pca_info.trim.enabled = do_trim; 
pca_info.trim.idx_kept = idx_kept;

pca_info.medfilt.enabled = do_medfilt;  %#ok<*STRNU>
pca_info.medfilt.halfwidth = medfilt_halfwidth;

save(savename, 'pca_info', 'pca_filters', 'pca_traces', 'S');

fprintf('%s: All done!\n', datestr(now));
function run_pca(movie_source, varargin)
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
autoset_num_PCs = 1;
max_num_PCs = 1000;


if ~isempty(varargin)
    for k = 1:length(varargin)
        switch varargin{k}
            case 'trim'
                do_trim = 1;
            case 'medfilt'
                do_medfilt = 1;
            case 'num_PCs'
                num_PCs = varargin{k+1};
                if ~isinteger(num_PCs)
                    error('num_PCs must be an integer.');
                end
                autoset_num_PCs = 0;
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

% Calculate max-projection (required for trimming and auto-setting #of PCs)
if do_trim || autoset_num_PCs
    max_proj = max(M,[],3);
end

% Detect local maxima in max-projection image to get an estimated maximum 
%number of cells in the movie
if autoset_num_PCs
    cents = local_maxima_2D(max_proj);
    num_PCs = size(cents,2);
    if num_PCs>max_num_PCs
        num_PCs = max_num_PCs;
    end        
    fprintf('%s: Extracting %d PCs...\n', datestr(now),num_PCs);
end

% Reshape movie into [space x time] matrix
num_pixels = height * width;
M = reshape(M, num_pixels, num_frames);

% Make each frame zero-mean in place
mean_M = mean(M,1);
M = bsxfun(@minus, M, mean_M);

idx_kept = 1:num_pixels;
if do_trim
    idx_kept = find(max_proj(:)>median(max_proj(:)));
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
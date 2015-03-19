function reconstruct_traces(movie_source, ica_dir, varargin)
% Reconstruct trace from the movie using the IC filters.
%
% inputs:
%
%   movie_source : Name of the movie file
%   ica_dir: Directory containing ICA results in a "ica_*.mat" file
%
% variable input arguments:
%
%   threshMov : Scale for the lower threshold for the movie brightness.( Any
%           value lower than this threshold will be set to 0). If threshMov 
%           is negative, then the lower threshold is calculated as 
%           |threshMov|*min(movie_matrix). If it is positive, then the lower 
%           threshold is calculated as threshMov * max(movie_matrix).
%           Default for threshMov is 0. Setting it to -1 retains all movie
%           pixels.
%   threshIC : Scale for the threshold for the ICA filter weights. (any
%           value in an ica filter lower than the threshold will be set to 
%           0.) The threshold for an IC is calculated as 
%           threshIC * max(ica_filter).
%           Default for threshMov is 0.3 . 
%   
% output:
%
%   Mat file containing reconstructed filters and traces
%
% Example usage : 
%   reconstruct_traces('c9m7d06.hdf5','ica001','threshIC',0.2,'threshMov',-1)
%
% Hakan Inan (Mar 15)
%

% Reconstruction parameters
threshMov = -1; % By default keep everything
threshIC = 0.3;

if ~isempty(varargin)
    len = length(varargin);
    for k = 1:len
        switch varargin{k}
            case 'threshMov'
                threshMov = varargin{k+1};
                if ~isnumeric(threshMov)
                    error('Movie brightness threshold must be numerical');
                end
            case 'threshIC'
                threshIC = varargin{k+1};
                if ~isnumeric(threshIC)
                    error('IC filter threshold must be numerical');
                end
        end
    end
end

% Load data
%------------------------------------------------------------
fprintf('%s: Loading movie...\n', datestr(now));
M = load_movie(movie_source);
[height, width, num_frames] = size(M);
M = reshape(M, height*width, num_frames);

if threshMov<0
    minMov = min(M(:));
    threshMov = -threshMov*minMov;
elseif threshMov>0
    maxMov = max(M(:));
    threshMov = -threshMov*maxMov;
end

% Load ICA
ica_filename = get_most_recent_file(ica_dir, 'ica_*.mat');
ica = load(ica_filename);
[height, width, num_ICs] = size(ica.filters);

% Compute the reconstruction
%------------------------------------------------------------
info.type = 'reconstruction';
info.movie_source = movie_source;
info.ica_source = ica_filename;

info.threshIC = threshIC;
info.threshMov = threshMov;

filters = zeros(height, width, 2*num_ICs, 'single'); % Preallocate
fprintf('%s: Thresholding IC filters...\n', datestr(now));

rec_filter_count = 0;
for idx_cell = 1:num_ICs
    ic_filter = ica.filters(:,:,idx_cell);
    [boundaries, ic_mask] = compute_ic_boundary(ic_filter, threshIC);
    num_boundaries = length(boundaries);
    if (num_boundaries == 1) % Only 1 ROI detected
        rec_filter_count = rec_filter_count + 1;
        filters(:,:,rec_filter_count) = ic_filter .* ic_mask;
    else % Multiple ROIs detected -- prompt user for split
        imagesc(ic_filter .* ic_mask);
        hold on;
        for k = 1:num_boundaries
            boundary = boundaries{k};
            plot(boundary(:,1), boundary(:,2), 'r', 'LineWidth', 2);
            text(max(boundary(:,1)), min(boundary(:,2)),...
                 num2str(k),...
                 'Color', 'w');
        end
        hold off;
        
        title(sprintf('%s -- IC %d', strrep(ica_filename, '_', '\_'), idx_cell));
        axis image;
        pause;
    end
end
info.num_pairs = rec_filter_count; %#ok<STRNU>
filters = filters(:,:,1:rec_filter_count);

for idx_cell = 1:rec_filter_count % Normalize
    filter = filters(:,:,idx_cell);
    filters(:,:,idx_cell) = filter / sum(filter(:));
end
    
traces = zeros(num_frames, rec_filter_count, 'single');
fprintf('%s: Reconstructing traces...\n', datestr(now));
for idx_cell = 1:rec_filter_count
    rec_filter = filters(:,:,idx_cell);
    pix_active = find(rec_filter>0);
    movie_portion = M(pix_active,:)';
    movie_portion(movie_portion<threshMov) = 0;
    traces(:,idx_cell) = movie_portion * rec_filter(pix_active);  
end

% Save the result to mat file
%------------------------------------------------------------
timestamp = datestr(now, 'yymmdd-HHMMSS');
rec_savename = sprintf('rec_%s.mat', timestamp);

save(rec_savename, 'info', 'filters', 'traces');

fprintf('%s: Done!\n', datestr(now));


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
%Hakan Inan (Mar 15)
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

% Compute the reconstruction
%------------------------------------------------------------
info.type = 'reconstruction';
info.movie_source = movie_source;
info.ica_source = ica_filename;

info.num_pairs = ica.info.num_ICs; % In the future, may be different
info.threshIC = threshIC;
info.threshMov = threshMov; %#ok<STRNU>

num_ICs = ica.info.num_ICs;
fprintf('%s: Thresholding IC filters...\n', datestr(now));
filters = zeros(size(ica.filters), 'single');
for idx_cell = 1:num_ICs
    filters(:,:,idx_cell) = threshold_ic_filter(ica.filters(:,:,idx_cell),threshIC);
end

traces = zeros(num_frames, num_ICs, 'single');
fprintf('%s: Reconstructing traces...\n', datestr(now));
for idx_cell = 1:num_ICs
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


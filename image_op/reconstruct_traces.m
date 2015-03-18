function reconst_traces = reconstruct_traces(movie_source,ica_filters,varargin)
% Reconstruct trace from the movie using the IC filters.
%
% inputs:
%
%   movie_source : link to the movie file
%   ica_filters : ica filters corresponding to the movie
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
%   reconst_traces = [# of frames] x [# if ICs] array of reconstructed
%   traces
%
% Example usage : 
% reconst_traces = reconstruct_traces(movie_source,ica_filters,'threshIC',0.2,'threshMov',-1)
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

fprintf('%s: Loading movie...\n', datestr(now));
M = load_movie(movie_source);
[height,width,num_frames] = size(M);
M = reshape(M,height*width,num_frames);

if threshMov<0
    minMov = min(M(:));
    threshMov = -threshMov*minMov;
elseif threshMov>0
    maxMov = max(M(:));
    threshMov = -threshMov*maxMov;
end

num_ICs = size(ica_filters,3);
fprintf('%s: Thresholding IC filters...\n', datestr(now));
for idx_cell = 1:num_ICs
    ica_filters(:,:,idx_cell) = threshold_ic_filter(ica_filters(:,:,idx_cell),threshIC);
end

reconst_traces = zeros(num_frames,num_ICs);
fprintf('%s: Reconstructing traces...\n', datestr(now));
for idx_cell = 1:num_ICs
    ica_filter = ica_filters(:,:,idx_cell);
    pix_active = find(ica_filter>0);
    movie_portion = M(pix_active,:)';
    movie_portion(movie_portion<threshMov) = 0;
    reconst_traces(:,idx_cell) = movie_portion * ica_filter(pix_active);  
end
fprintf('%s: Done!\n', datestr(now));


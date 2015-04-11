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

% Minimum number of active pixels in an IC filter to be regarded as valid
min_num_pixels = 4;

% Reconstruction parameters
threshMov = -1; 
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

% Load ICA
ica_filename = get_most_recent_file(ica_dir, 'ica_*.mat');
ica = load(ica_filename);
[height, width, num_ICs] = size(ica.filters);

% Reconstruction settings
info.type = 'reconstruction';
info.movie_source = movie_source;
info.ica_source = ica_filename;

info.threshIC = threshIC;
info.threshMov = threshMov;

% Build reconstruction filters
%------------------------------------------------------------
filters = zeros(height, width, 2*num_ICs, 'single'); % Preallocate buffer
fprintf('%s: Thresholding IC filters...\n', datestr(now));

rec_filter_count = 0;
for ic_idx = 1:num_ICs
    ic_filter = ica.filters(:,:,ic_idx);
    [boundaries, ic_mask] = compute_ic_boundary(ic_filter, threshIC);
    num_boundaries = length(boundaries);
    if (num_boundaries == 1) % Only 1 ROI detected
        rec_filter_count = rec_filter_count + 1;
        filters(:,:,rec_filter_count) = ic_filter .* ic_mask;
    else % Multiple ROIs detected
        if (num_boundaries > 20)
            fprintf('%s: IC %d has %d boundaries, skipping...\n',...
                datestr(now), ic_idx, num_boundaries);
            continue;
        end
        
        submasks = zeros(height, width, num_boundaries);
        
        % Show the multiple ROIs to the user
        imagesc(ic_filter .* ic_mask);
        colormap gray;
        hold on;
        for k = 1:num_boundaries
            if mod(k,2)
                color = 'g';
            else
                color = 'r';
            end
            boundary = boundaries{k};
            plot(boundary(:,1), boundary(:,2), color, 'LineWidth', 1);
            text(max(boundary(:,1)), min(boundary(:,2)),...
                 num2str(k),...
                 'Color', color);
             
            submasks(:,:,k) = poly2mask(boundary(:,1), boundary(:,2), height, width);
        end
        hold off;
        
        title(sprintf('%s -- IC %d', ica_filename, ic_idx),...
              'Interpreter', 'none'); % Turn off Tex
        axis image;
        drawnow;
        
        % Get list of ROIs to keep
        fprintf('%s: Please select ROIs to keep from IC %d\n', datestr(now), ic_idx);
        sel_rois = [];
        val = str2double(strtrim(input('  >> ','s')));
        while (1)
            if ~isnan(val) % Is a number
                if ((1 <= val) && (val <= num_boundaries))
                    if ismember(val, sel_rois)
                        sel_rois = setdiff(sel_rois, val);
                    else
                        sel_rois = union(sel_rois, val);
                    end
                else
                    fprintf('    Sorry, %d is not a valid ROI index\n', val);
                end
            else % Not a number -- done with ROI selection
                for sel_roi = sel_rois
                    rec_filter_count = rec_filter_count + 1;
                    filters(:,:,rec_filter_count) = ic_filter .* submasks(:,:,sel_roi);
                    fprintf('    Added ROI %d from IC %d\n', sel_roi, ic_idx);
                end
                break; % Get out of selection loop
            end
            
            val = str2double(strtrim(input('  >> ','s')));
        end
    end
end
filters = filters(:,:,1:rec_filter_count);

cells_to_exclude = [];
for cell_idx = 1:rec_filter_count % Normalize
    filter = filters(:,:,cell_idx);
    sum_filter = sum(filter(:));
    if sum_filter>0 && sum(filter(:)>0)>=min_num_pixels
        filters(:,:,cell_idx) = filter / sum_filter;
    else
        cells_to_exclude(end+1) = cell_idx; 
    end
end
filters(:,:,cells_to_exclude) = []; 
rec_filter_count = rec_filter_count - length(cells_to_exclude);
info.num_pairs = rec_filter_count;%#ok<STRNU>

% Reconstruct traces
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

traces = zeros(num_frames, rec_filter_count, 'single');
fprintf('%s: Reconstructing traces...\n', datestr(now));
for cell_idx = 1:rec_filter_count
    rec_filter = filters(:,:,cell_idx);
    pix_active = find(rec_filter>0);
    movie_portion = M(pix_active,:)';
    movie_portion(movie_portion<threshMov) = 0;
    trace_this = (movie_portion * rec_filter(pix_active))';  
    traces(:,cell_idx) = fix_baseline(trace_this);
end

% Save the result to mat file
%------------------------------------------------------------
timestamp = datestr(now, 'yymmdd-HHMMSS');
rec_savename = sprintf('rec_%s.mat', timestamp);

save(rec_savename, 'info', 'filters', 'traces');

fprintf('%s: Done!\n', datestr(now));

function [rec_savename, class_savename] = backapply_filters(filters_in, movie_in, varargin)
% Compute traces by back applying filters to the specified movie. Note that
% output filters are normalized such that the sum of all pixels equals 1.
%
% Inputs:
%   - 'filters_in': Can be a DaySummary instance, or a 3D matrix of filters
%       where filters_in(:,:,k) is the image of the k-th filter.
%   - 'movie_in': Can be the name of a HDF5 file (string) or a 3D matrix.
%

use_ls = false;
use_fix_baseline = false;

generate_class = false;
class_savename = [];

% 'merge_list' allows multiple filters in 'filters_in' to be combined into
% a single filter, and should be formatted as:
%   merge_list = {[a b c], [d e], [f g h]}
% which specifies input filter #'s a/b/c should be combined into one 
% output filter, d/e into one output filter, etc.
merge_list = {};

for k = 1:length(varargin)
    vararg = varargin{k};
    if ischar(vararg)
        switch lower(vararg)
            case {'ls', 'leastsquares'}
                use_ls = true;
            case {'fix', 'fix_baseline'}
                use_fix_baseline = true;
            case {'class', 'generate_class'}
                generate_class = true;
            case 'merge'
                merge_list = varargin{k+1};
        end
    end
end

% Movie input can be a filename (HDF5) or a 3-D movie matrix
%------------------------------------------------------------
if ischar(movie_in)
    movie_dims = get_movie_info(movie_in);

    movie_dataset = '/Data/Images';
    get_frames = @(x) h5read(movie_in, movie_dataset, [1 1 x(1)], [movie_dims(1) movie_dims(2) x(end)-x(1)+1]);
else
    movie_dims = size(movie_in);
    
    get_frames = @(x) movie_in(:,:,x);
end
height = movie_dims(1);
width = movie_dims(2);
num_pixels = height * width;
num_frames = movie_dims(3);

% Build up 'filters'
%------------------------------------------------------------
fprintf('%s: Building filters matrix...\n', datestr(now));
if isa(filters_in, 'DaySummary')
    cell_indices = find(filters_in.is_cell);
    get_filter = @(x) filters_in.cells(x).im;
else
    cell_indices = 1:size(filters_in,3);
    get_filter = @(x) filters_in(:,:,x);
end

% Filters to be merged will be handled separately
num_merged_filters = length(merge_list);
filters_to_be_merged = cell2mat(merge_list);
cell_indices = setdiff(cell_indices, filters_to_be_merged);

num_filters = length(cell_indices); % Number of filters, excluding those to be merged
filters = zeros(num_pixels, num_filters + num_merged_filters, 'single'); % Preallocate

idx = 0;
num_skipped = 0;
for k = 1:num_filters
    cell_idx = cell_indices(k);
    filter = get_filter(cell_idx);

    filter_sum = sum(filter(:)); % Would it be better to take the absolute value?
    if (filter_sum > 0)
        idx = idx + 1;
        filter = filter / filter_sum; % Normalize
        filters(:,idx) = reshape(filter, num_pixels, 1);
    else
        fprintf('  Warning: Filter #%d is entirely zero -- skipping!\n', k);
        num_skipped = num_skipped + 1;
    end
end

% Generate merged filters
for k = 1:num_merged_filters
    merge_inds = merge_list{k};
    merge_inds3 = permute(merge_inds, [1 3 2]); % Inds stacked along 3rd dim of array
    fs = arrayfun(get_filter, merge_inds3, 'UniformOutput', false);
    fs = cell2mat(fs); % Filters are stacked along 3rd dim of 'fs'
    
    filter = sum(fs,3);
    filter_sum = sum(filter(:));
    if (filter_sum > 0)
        idx = idx + 1;
        filter = filter / filter_sum;
        filters(:,idx) = reshape(filter, num_pixels, 1);
        fprintf('  Merged filters [%s] into one ouput filter\n', num2str(merge_inds));
    else
        fprintf('  Warning: Merge of filters [%s] is entirely zero -- skipping!\n', num2str(merge_inds));
        num_skipped = num_skipped + 1;
    end
end

if (num_skipped > 0)
    fprintf('  Total number of skipped filters: %d\n', num_skipped);
end

filters(:,(idx+1):end) = [];
num_filters = idx;

% Calculate traces in chunks
%------------------------------------------------------------
traces = zeros(num_filters, num_frames, 'single');

frame_chunk_size = 1000;
[frame_chunks, num_chunks] = make_frame_chunks(num_frames, frame_chunk_size);

for i = 1:num_chunks
    frame_inds = frame_chunks(i,1):frame_chunks(i,2); % For this chunk
    num_frames_in_chunk = frame_inds(end) - frame_inds(1) + 1;
    
    fprintf('%s: Calculating traces for frames %d to %d (out of %d) ',...
        datestr(now), frame_inds(1), frame_inds(end), num_frames);
        
    movie_chunk = get_frames(frame_inds);
    movie_chunk = reshape(movie_chunk, num_pixels, num_frames_in_chunk);
    
    %------------------------------------------------------------
    % Comments regarding trace scaling (2021 Aug 15):
    %
    %   Given the probability-like normalization of filter (sum to 1)
    %   above, "simple projection" is effectively taking a weighted average
    %   across pixels of the movie. Thus, the projection trace will have
    %   units and numerical scaling that is directly comparable to the
    %   movie itself.
    %
    %   Note, however, that when filters are normalized as probabilities,
    %   the typical values in each pixel will be << 1, since usually many
    %   pixels will be nonzero in the filter, and the sum must equal 1.
    %
    %   This filter normalization causes numerical values in the least
    %   squares trace to be a lot larger than the analogous simple
    %   projection trace. This is because the least squares fit is trying
    %   to reconstruct the movie frame using the provided filter whose
    %   values are << 1. Thus the trace values need to "compensate" the
    %   small numerical values of the spatial filter.
    %
    %   To _roughly_ re-scale the least square trace to the scaling of the
    %   input movie, a procedure like the following can be used. Let:
    %       f: 2D spatial filter of k-th cell, linearized as a vector
    %       t: least-squares trace
    %   then,
    %       fm = max(f(:));
    %       t2 = mean(f(f>0.3*fm))) * t;
    %   Here, t2 will have scaling comparable to the input movie.
    %
    %   The basic idea is that, for least-squares to produce trace values
    %   comparable to the original movie, we need to normalize the
    %   _average_ value of the "active" pixels in the filter to be 1,
    %   rather than the sum of the filter pixels. The threshold f > 0.3*fm
    %   is a heuristic for selecting the "active" pixels in the filter.
    %------------------------------------------------------------
    if use_ls
        fprintf('with least squares...\n');
        traces(:,frame_inds) = filters \ movie_chunk;
    else % Otherwise, simple projection
        fprintf('with simple projection...\n');
        traces(:,frame_inds) = filters' * movie_chunk;
    end
end

% Reshape to standard form
filters = reshape(filters, height, width, num_filters);
traces = traces'; % [num_frames x num_filters]

% Post-process traces
%------------------------------------------------------------
if use_fix_baseline
    fprintf('%s: Applying trace baseline (running percentile) correction... ', datestr(now));
    tic;
    for k = 1:num_filters
        traces(:,k) = fix_baseline(traces(:,k));
    end
    t = toc;
    fprintf('Done (%.1f s)\n', t);
end

% Save as Rec file
%------------------------------------------------------------
if ischar(movie_in)
    info.type = sprintf('backapply_filters:%s', movie_in);
else
    info.type = 'backapply_filters';
end
info.num_pairs = num_filters;

% Note the parameters used in DFF recomputation
info.options.use_ls = use_ls;
info.options.use_fix_baseline = use_fix_baseline;

tic;
[rec_savename, timestamp] = save_rec(info, filters, traces);
if generate_class
    class_savename = generate_class_file(num_filters, 'timestamp', timestamp);
end
t = toc;
fprintf('%s: %d filter-trace pairs saved to "%s" (%.1f sec)\n', datestr(now), num_filters, rec_savename, t);
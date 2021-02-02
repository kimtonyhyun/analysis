function rec_savename = backapply_filters(ds, movie_in, varargin)
% Compute traces by back applying filters from DaySummary to the specified
% movie file. Originally based on 'get_dff_traces'.

use_ls = false;

% Filter parameters relevant for trace extraction
use_all_filters = false;
truncate_filter = false;

% Trace post-extraction processing
fix_baseline_method = [];

for k = 1:length(varargin)
    vararg = varargin{k};
    if ischar(vararg)
        switch lower(vararg)
            case 'keepall'
                use_all_filters = true;
            case {'ls', 'leastsquares'}
                use_ls = true;
            case {'fix', 'fix_baseline'}
                fix_baseline_method = varargin{k+1};
            case 'truncate'
                fprintf('%s: Filter will be truncated...\n', datestr(now));
                truncate_filter = true;
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
fprintf('%s: Building filters matrix... ', datestr(now));
tic;
if use_all_filters
    cell_indices = 1:ds.num_cells;
else 
    cell_indices = find(ds.is_cell);
end
num_filters = length(cell_indices);
filters = zeros(num_pixels, num_filters, 'single');

for k = 1:num_filters
    cell_idx = cell_indices(k);
    filter = ds.cells(cell_idx).im;
    if truncate_filter
        filter = filter .* ds.cells(cell_idx).mask;
    end
    % NOTE: In principle, we should check if the filter is entirely 0
    % here, as in the section below.
    filter = filter / sum(filter(:)); % Normalize filter to 1

    filters(:,k) = reshape(filter, num_pixels, 1);
end
t = toc;
fprintf('Done (%.1f sec)\n', t);

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
    
    if use_ls
        fprintf('with least squares...\n');
        traces(:,frame_inds) = filters \ movie_chunk;
    else % Otherwise, simple projection
        fprintf('with simple projection...\n');
        traces(:,frame_inds) = filters' * movie_chunk;
    end
end

% Post-process traces
%------------------------------------------------------------
if ~isempty(fix_baseline_method)
    fprintf('%s: Applying baseline correction (%s) to DFF traces... ', datestr(now), fix_baseline_method);
    tic;
    for k = 1:num_filters
        traces(k,:) = fix_baseline(traces(k,:), fix_baseline_method);
    end
    t = toc;
    fprintf('(%.1f s)\n', t);
end

% Save as Rec file
%------------------------------------------------------------

% Reshape to standard form
filters = reshape(filters, height, width, num_filters);
traces = traces'; % [num_frames x num_filters]

% Save as Rec file
%------------------------------------------------------------
if ischar(movie_in)
    info.type = sprintf('backapply_filters:%s', movie_in);
else
    info.type = 'backapply_filters';
end
info.num_pairs = num_filters;

% Note the parameters used in DFF recomputation
info.options.fix_baseline = fix_baseline_method;
info.options.truncate_filter = truncate_filter;
info.options.use_ls = use_ls;
info.options.use_all_filters = use_all_filters;

tic;
rec_savename = save_rec(info, filters, traces);
t = toc;
fprintf('%s: Results saved to "%s" (%.1f sec)\n', datestr(now), rec_savename, t);
function [traces, rec_savename] = get_dff_traces(filter_input, M_dff, varargin)
% Compute DFF traces associated with classified cells from a DaySummary,
% and saves the result as a rec file. Assumes that M_dff is a DFF movie (in
% order for the resulting traces to be interpreted as DFF).
%
% Note that the length of traces is derived from the length of the movie.
% The DaySummary only provides the spatial filters.
%
% TODO:
%   - Load chunks of the movie at a time to reduce memory footprint

use_ls = false;
suppress_output = false;

% Filter parameters relevant for trace extraction -- used only when
% 'filter_input' is a DaySummary instance
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
            case 'no_rec'
                fprintf('get_dff_traces: No output file will be generated.\n');
                suppress_output = true;
        end
    end
end

[height, width, num_frames] = size(M_dff);
num_pixels = height*width;
M_dff = reshape(M_dff, num_pixels, num_frames); % Reshape movie as a matrix

if isa(filter_input, 'DaySummary')
    % Filter input is provided as a DaySummary object
    ds = filter_input;

    if use_all_filters
        % Use all filters regardless of cell classification
        cell_indices = 1:ds.num_cells;
    else 
        % Use only filters that have been classified to be a cell
        cell_indices = find(ds.is_cell);
    end
    num_filters = length(cell_indices);
    filters = zeros(num_pixels, num_filters);

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
else
    % Otherwise, we assume 'filter_input' is a [height x width x
    % num_filters] matrix of filter images.
    num_filters = size(filter_input, 3);
    filters = zeros(num_pixels, num_filters);
    
    idx = 0;
    for k = 1:num_filters
        filter = filter_input(:,:,k);
        if (sum(filter(:)) > 0)
            filter = filter / sum(filter(:));
            
            idx = idx + 1;
            filters(:,idx) = reshape(filter, num_pixels, 1);
        else
            fprintf('Warning: Filter #%d is entirely zero -- skipping!\n', k);
        end
    end

    num_filters = idx;
    filters = filters(:,1:idx);    
end
fprintf('%s: Computing traces ', datestr(now));
tic;
if use_ls
    fprintf('by least squares... ');
    traces = filters \ M_dff;
    traces = double(traces); % 'traces' takes the type of 'M_dff', which is typically single
else % Simple projection
    fprintf('by simple projection... ');
    traces = filters' * M_dff;
end
t = toc;
fprintf('Done! (%.1f s)\n', t);

% Reshape to standard form
filters = reshape(filters, height, width, num_filters);
traces = traces'; % [num_frames x num_cells]

if ~isempty(fix_baseline_method)
    fprintf('%s: Applying baseline correction (%s) to DFF traces... ', datestr(now), fix_baseline_method);
    tic;
    for k = 1:num_filters
        trace = traces(:,k);
        traces(:,k) = fix_baseline(trace, fix_baseline_method);
    end
    t = toc;
    fprintf('(%.1f s)', t);
end

% Save as Rec file
%------------------------------------------------------------
info.type = 'get_dff_traces';
info.num_pairs = num_filters;

% Note the parameters used in DFF recomputation
info.options.fix_baseline = fix_baseline_method;
info.options.truncate_filter = truncate_filter;
info.options.use_ls = use_ls;
info.options.use_all_filters = use_all_filters;

rec_savename = [];
if ~suppress_output
    rec_savename = save_rec(info, filters, traces);
end
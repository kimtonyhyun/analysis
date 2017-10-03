function traces = get_dff_traces(ds, M_dff, varargin)
% Compute DFF traces associated with classified cells from a DaySummary,
% and saves the result as a rec file. Assumes that M_dff is a DFF movie!
%
% Note that the length of traces is derived from the length of the movie.
% The DaySummary only provides the spatial filters.

truncate_filter = false;
remove_baseline = false;
use_ls = false;

for k = 1:length(varargin)
    vararg = varargin{k};
    if ischar(vararg)
        switch lower(vararg)
            case {'ls', 'leastsquares'}
                use_ls = true;
            case 'fix_baseline'
                remove_baseline = true;
            case 'truncate'
                fprintf('%s: Filter will be truncated...\n', datestr(now));
                truncate_filter = true;
        end
    end
end

[height, width, num_frames] = size(M_dff);
num_pixels = height*width;
M_dff = reshape(M_dff, num_pixels, num_frames); % Reshape movie as a matrix

cell_indices = find(ds.is_cell);
num_cells = length(cell_indices);
filters = zeros(num_pixels, num_cells); % Filter of actual cells only

fprintf('%s: Building filter matrix of %d classified cells...\n',...
    datestr(now), num_cells);
for k = 1:num_cells
    cell_idx = cell_indices(k);
    filter = ds.cells(cell_idx).im;
    if truncate_filter
        filter = filter .* ds.cells(cell_idx).mask;
    end
    filter = filter / sum(filter(:)); % Normalize filter to 1
    
    filters(:,k) = reshape(filter, num_pixels, 1);
end

fprintf('%s: Computing traces... ', datestr(now));
tic;
if use_ls
    traces = filters \ M_dff;
else % Simple projection
    traces = filters' * M_dff;
end
t = toc;
fprintf('Done! (%.1f s)\n', t);

% Reshape to standard form
filters = reshape(filters, height, width, num_cells);
traces = traces'; % [num_frames x num_cells]

if remove_baseline
    fprintf('%s: Applying baseline correction to DFF traces...\n', datestr(now));
    for k = 1:num_cells
        trace = traces(:,k);
        traces(:,k) = fix_baseline(trace);
    end
end

% Save as Rec file
%------------------------------------------------------------
info.type = 'get_dff_traces';
info.num_pairs = num_cells;

% Note the parameters used in DFF recomputation
info.options.remove_baseline = remove_baseline;
info.options.truncate_filter = truncate_filter;

timestamp = datestr(now, 'yymmdd-HHMMSS');
rec_savename = sprintf('rec_%s.mat', timestamp);
save(rec_savename, 'info', 'filters', 'traces', '-v7.3');
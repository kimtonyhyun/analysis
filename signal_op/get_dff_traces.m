function traces = get_dff_traces(ds, M_dff)
% Compute DFF traces associated with classified cells from a DaySummary,
% and saves the result as a rec file. Assumes that M_dff is a DFF movie!

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
    filter = filter / sum(filter(:)); % Normalize filter to 1
    
    filters(:,k) = reshape(filter, num_pixels, 1);
end

fprintf('%s: Computing traces...\n', datestr(now));
traces = filters' * M_dff;

% Save as Rec file
%------------------------------------------------------------
info.type = 'get_dff_traces';
info.num_pairs = num_cells;

filters = reshape(filters, height, width, num_cells);
traces = traces'; % [num_frames x num_pairs]

timestamp = datestr(now, 'yymmdd-HHMMSS');
rec_savename = sprintf('rec_%s.mat', timestamp);
save(rec_savename, 'info', 'filters', 'traces', '-v7.3');
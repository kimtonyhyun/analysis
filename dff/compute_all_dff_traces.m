function [rec_savename, class_savename] = compute_all_dff_traces(ds, varargin)
% Compute DFF traces for all cells in DaySummary. Applies the same baseline
% estimation parameters to all cells.

num_frames = ds.full_num_frames;
[height, width] = size(ds.cell_map_ref_img);

cell_indices = find(ds.is_cell);
num_cells = length(cell_indices);

filters = zeros(height, width, num_cells, 'single');
traces = zeros(num_frames, num_cells, 'single');
dff_infos = cell(num_cells, 1);

for k = 1:num_cells
    if mod(k,25) == 1
        fprintf('%s: Computing DFF traces for cell %d of %d...\n', datestr(now), k, num_cells);
    end
    cell_idx = cell_indices(k);
    filters(:,:,k) = ds.cells(cell_idx).im;
    
    tr = ds.get_trace(cell_idx);
    [traces(:,k), dff_infos{k}] = compute_dff_trace(tr, varargin{:});
end
fprintf('%s: Done!\n', datestr(now));

% Save as Rec file
%------------------------------------------------------------
info.type = 'compute_all_dff_traces';
info.num_pairs = num_cells;

info.dff_infos = cell2mat(dff_infos);

[rec_savename, timestamp] = save_rec(info, filters, traces);
class_savename = generate_class_file(num_cells, 'timestamp', timestamp);
fprintf('%s: %d filter-trace pairs saved to "%s"\n', datestr(now), num_cells, rec_savename);
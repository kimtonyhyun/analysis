function [rec_savename, class_savename] = get_dff_traces(ds_proj, ds_ls, fps)
% TODO:
%   - Build out a save/load mechanism
%

num_cells = ds_ls.num_cells;
num_frames = ds_ls.full_num_frames;

state.fig_handle = figure;
state.baseline_method = 2;
state.running_percentile_window = 100 * fps; % Frames

% Populate the following variables
dff_traces = zeros(num_frames, num_cells, 'single');
dff_infos = cell(num_cells, 1);

% Flag where:
%   - keep_dff{k} == true;  % Keep the k-th cell in export
%   - keep_dff{k} == false; % Omit the k-th cell in export
%   - keep_dff{k} == [];    % K-th cell not yet classified
keep_dff = cell(num_cells, 1);

% Application loop
%------------------------------------------------------------
cell_idx = 1;
prev_cell_idx = 1;
while (cell_idx <= num_cells)
    switch state.baseline_method
        case 0 % Mode
            [dff_trace, dff_info] = inspect_dff_traces(ds_proj, ds_ls, cell_idx, fps, 'mode');
        case 1 % Linear
            [dff_trace, dff_info] = inspect_dff_traces(ds_proj, ds_ls, cell_idx, fps, 'linear');
        case 2 % Running percentile
            [dff_trace, dff_info] = inspect_dff_traces(ds_proj, ds_ls, cell_idx, fps,...
                'percentile', 'window', state.running_percentile_window);
    end
    
    % Ask the user to classify the cell candidate
    prompt = sprintf('Get DFF traces (%d/%d; bmode=%d) >> ', cell_idx, num_cells, state.baseline_method);
    resp = strtrim(input(prompt, 's'));
    
    val = str2double(resp);
    if (~isnan(val)) % Is a number. Check if it is a valid index and jump to it
        if ((1 <= val) && (val <= num_cells))
            prev_cell_idx = cell_idx;
            cell_idx = val;
        else
            fprintf('  Sorry, %d is not a valid cell index\n', val);
        end
    else
        resp = lower(resp);
        switch (resp(1))
            case 'c' % Accept DFF result
                dff_traces(:,cell_idx) = dff_trace;
                dff_infos{cell_idx} = dff_info;
                keep_dff{cell_idx} = true;
                
                go_to_next_cell();
                
            case 'n' % Reject DFF
                fprintf('  Cell %d will be omitted from DFF output\n', cell_idx);
                keep_dff{cell_idx} = false;
                go_to_next_cell();
                
            case ''
                go_to_next_cell();
                
            case 'b' % Select baseline correction mode
                switch resp
                    case 'b0'
                        state.baseline_method = 0;
                        fprintf('  Baseline is computed by the trace mode\n');
                    case 'b1'
                        state.baseline_method = 1;
                        fprintf('  Baseline is computed via a linear fit\n');
                    case 'b2'
                        state.baseline_method = 2;
                        fprintf('  Baseline is computed as a running percentile\n');
                    otherwise
                        cprintf('red', '  Invalid baseline command\n');
                end

            case 'p' % Jump to previously viewed cell
                temp_idx = cell_idx;
                cell_idx = prev_cell_idx;
                prev_cell_idx = temp_idx;

            case 'q' % Exit
                close(state.fig_handle);
                break;
        end
    end
end

% Save as Rec file
%------------------------------------------------------------

% First, collect ALL spatial filters and traces
[height, width] = size(ds_ls.cell_map_ref_img);
filters = zeros(height, width, num_cells, 'single');
for k = 1:num_cells
    filters(:,:,k) = ds_ls.cells(k).im;
end

% Unlabeled entries in 'keep_dff' will NOT be exported
for k = 1:num_cells
    if isempty(keep_dff{k})
        keep_dff{k} = false;
    end
end
keep_dff = cell2mat(keep_dff);
num_exported_cells = sum(keep_dff);

filters = filters(:,:,keep_dff);
traces = dff_traces(:,keep_dff);

info.type = 'get_dff_traces';
info.num_pairs = num_exported_cells;

% Note that 'dff_infos' must remain a (Matlab) cell, since the contents of 
% the struct may be different across different cells in the DaySummary.
info.dff_infos = dff_infos;
info.fps = fps;

[rec_savename, timestamp] = save_rec(info, filters, traces);
class_savename = generate_class_file(num_exported_cells, 'timestamp', timestamp);
fprintf('%s: %d filter-trace pairs saved to "%s"\n', datestr(now), num_exported_cells, rec_savename);

    % Auxiliary functions
    %------------------------------------------------------------
    function go_to_next_cell()
        unprocessed = cellfun(@isempty, keep_dff);
        next_idx = find_next_cell_to_process(cell_idx, unprocessed);
        
        if isempty(next_idx)
            cprintf([0 0.5 0], '  All DFF traces have been computed!\n');
            prev_cell_idx = cell_idx;
            cell_idx = cell_idx + 1;
            if cell_idx > num_cells
                cell_idx = 1;
            end
        else
            prev_cell_idx = cell_idx;
            cell_idx = next_idx;
        end
        
        state.baseline_method = 2;
    end % go_to_next_cell
    
end % get_dff_traces
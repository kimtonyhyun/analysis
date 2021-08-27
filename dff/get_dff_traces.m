function [rec_savename, class_savename] = get_dff_traces(ds_proj, ds_ls, fps)

num_cells = ds_proj.num_cells;

state.fig_handle = figure;
state.fix_wandering_baseline = true;
state.running_percentile_window = 100 * fps; % Frames

% Populate the following variables
dff_traces = cell(num_cells, 1);
dff_infos = cell(num_cells, 1);

% Application loop
%------------------------------------------------------------
cell_idx = 1;
prev_cell_idx = 1;
while (cell_idx <= num_cells)
    if state.fix_wandering_baseline
        % Note: 'inspect_dff_traces' returns the DFF trace as computed from
        %   the least squares trace.
        [dff_trace, dff_info] = inspect_dff_traces(ds_proj, ds_ls, cell_idx, fps,...
            'fix_baseline', state.running_percentile_window);
    else
        [dff_trace, dff_info] = inspect_dff_traces(ds_proj, ds_ls, cell_idx, fps);
    end
    
    % Ask the user to classify the cell candidate
    prompt = sprintf('Get DFF traces (%d/%d) >> ', cell_idx, num_cells);
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
        switch (resp)
            case 'c' % Accept DFF result
                dff_traces{cell_idx} = dff_trace;
                dff_infos{cell_idx} = dff_info;
                
                state.fix_wandering_baseline = true;
                go_to_next_cell();
            case ''
                go_to_next_cell();
            case 'b' % Toggle wandering baseline
                state.fix_wandering_baseline = ~state.fix_wandering_baseline;
                if (state.fix_wandering_baseline)
                    fprintf('  Computing DFF trace with respect to a WANDERING baseline\n');
                else
                    fprintf('  Computing DFF trace with respect to a FIXED baseline\n');
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
% Extract spatial filters
[height, width] = size(ds_ls.cell_map_ref_img);
filters = zeros(height, width, num_cells, 'single');
for k = 1:num_cells
    filters(:,:,k) = ds_ls.cells(k).im;
end
traces = cell2mat(dff_traces)';

info.type = 'get_dff_traces';
info.num_pairs = num_cells;

info.dff_infos = dff_infos;

[rec_savename, timestamp] = save_rec(info, filters, traces);
class_savename = generate_class_file(num_cells, 'timestamp', timestamp);
fprintf('%s: %d filter-trace pairs saved to "%s"\n', datestr(now), num_cells, rec_savename);

    % Auxiliary functions
    %------------------------------------------------------------
    function go_to_next_cell()
        unprocessed = cellfun(@isempty, dff_traces);
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
    end % go_to_next_cell
    
end % get_dff_traces
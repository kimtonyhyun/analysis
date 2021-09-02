function [rec_savename, class_savename] = get_dff_traces(ds_proj, ds_ls, fps)
% TODO:
%   - Build out a save/load mechanism
%

num_cells = ds_ls.num_cells;
num_frames = ds_ls.full_num_frames;

state.fig_handle = figure;

state.nonactive_thresh = 0.5;
state.nonactive_padding = 10 * fps;
state.polyfit_order = 1;

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
while (cell_idx <= num_cells)

    [dff_trace, dff_info, ax2] = inspect_dff_traces(ds_proj, ds_ls, cell_idx, fps,...
        'nonactive_polyfit', 'order', state.polyfit_order,...
        'thresh', state.nonactive_thresh,...
        'padding', state.nonactive_padding);
    
    % Ask the user to classify the cell candidate
    prompt = sprintf('Get DFF traces (%d/%d; order=%d) >> ', cell_idx, num_cells, state.polyfit_order);
    resp = strtrim(input(prompt, 's'));
    
    val = str2double(resp);
    if (~isnan(val)) % Is a number. Check if it is a valid index and jump to it
        if ((1 <= val) && (val <= num_cells))
            cell_idx = val;
        else
            fprintf('  Sorry, %d is not a valid cell index\n', val);
        end
    else
        resp = lower(resp);
        if isempty(resp)
            go_to_next_cell();
        else
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

                % Set baseline fitting parameters
                %------------------------------------------------------------                               
                case 'o' % Set order
                    try
                        val = str2double(resp(2:end));
                        if (~isnan(val)) % Is a number
                            if (val > 0)
                                state.polyfit_order = val;
                            end
                        end
                    catch
                        cprintf('red', '  Could not parse order command!\n');
                    end

                case 't' % Set threshold
                    fprintf('  Plese select a new threshold on the LS trace\n');
                    while (1)
                        [~, thresh] = ginput(1);
                        if (gca == ax2)
                            break;
                        else
                            fprintf('  Error! New threshold must be defined on the LS trace\n');
                        end
                    end

                    tr = ds_ls.get_trace(cell_idx);
                    state.nonactive_thresh = (thresh - min(tr))/(max(tr) - min(tr));
                    fprintf('  New activity threshold is %.1f%% of max\n', state.nonactive_thresh * 100);

                case 'q' % Exit
                    close(state.fig_handle);
                    break;
            end % switch(resp)
        end % isempty(resp)
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
info.dff_infos = dff_infos(keep_dff);
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
            cell_idx = cell_idx + 1;
            if cell_idx > num_cells
                cell_idx = 1;
            end
        else
            cell_idx = next_idx;
        end
        
        % Reset default fitting method on advance
        state.baseline_method = 'polyfit';
        state.nonactive_thresh = 0.5;
        state.polyfit_order = 1;
    end % go_to_next_cell
    
end % get_dff_traces
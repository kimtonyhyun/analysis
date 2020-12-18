function remove_duplicates(ds1, ds2)
% Given cell maps in two planes, look for duplicates by considering
% spatiotemporal alignment of cells, and then "assign" the cell to one of
% the DaySummary's. Note that the labels are mutated _in place_.

app_name = 'Remove duplicates';

% Compute all temporal correlations. Results are sorted by the correlation
% value (descending).
corrlist = compute_corrlist(ds1, ds2);
num_cells1 = ds1.num_classified_cells;
num_cells2 = ds2.num_classified_cells;
num_pairs = size(corrlist, 1);

% Sort corrlist by ds1 cell index. For a given ds1 cell, the corrlist rows
% then enumerate ds2 cells in order of decreasing correlation.
corrlist = sortrows(corrlist, 1, 'ascend');
idx1_to_i = corrlist(1:num_cells2:end,1);

hf = figure;

% Note that 'idxX' are not the actual cell indices into dsX. These are
% counters for traversing the corrlist.
idx1 = 1;
idx2 = 1;

while (1)
    % i and j are the actual cell indices into ds1 and ds2
    i = idx1_to_i(idx1);
    
    % Block of corrlist for the i-th cell in ds1
    rows = (1+(idx1-1)*num_cells2):(idx1*num_cells2);
    corrlist_i = corrlist(rows,:);
    
    j = corrlist_i(idx2,2);
    c = corrlist_i(idx2,3);
    show_corr(ds1, i, ds2, j, c, 'zsc', 'overlay', [], 'zoom_target', 1);
    
    prompt = sprintf('%s (ds1 idx1=%d of %d; ds2 idx2=%d of %d) >> ',...
                      app_name,...
                      idx1, num_cells1,...
                      idx2, num_cells2);
    resp = strtrim(input(prompt, 's'));
    
    val = str2double(resp);
    if (~isnan(val)) % Is a number
        if (1 <= val) && (val <= num_pairs)
            set_idx1(val);
        end
    else
        resp = lower(resp);
        if isempty(resp)
            % Increment ds2 cell index
            idx2 = idx2 + 1;
            idx2 = min(idx2, num_cells2);
        else
            switch resp(1)
                case 'c' % Jump to a cell in ds1 using its cell index
                    val = str2double(resp(2:end));
                    new_idx1 = find(idx1_to_i == val, 1, 'first');
                    if ~isempty(new_idx1)
                        set_idx1(new_idx1);
                    end
                    
                case 'a' % "Assign"
                    val = str2double(resp(2:end));
                    switch val
                        case 1 % Assign to ds1
                            ds2.cells(j).label = 'not a cell';
                            fprintf('  Assigned to ds1!\n');
                            increment_i;
                            
                        case 2 % Assign to ds2
                            ds1.cells(i).label = 'not a cell';
                            fprintf('  Assigned to ds2!\n');
                            increment_i;
                            
                        otherwise
                            cprintf('red', '  Error: Valid inputs are "a1" or "a2"\n');
                    end

                case 'm' % Cell map
                    display_map(i);

                case 'n' % Next
                    increment_i;
                    
                case 'p' % Previous
                    decrement_i;

                case 'q' % Exit
                    close(hf);
                    break;

                otherwise
                    fprintf('  Could not parse "%s"\n', resp);
            end
        end
    end
end

    % Display functions
    %------------------------------------------------------------
    function display_map(cell_idx)
        clf;
        color_mappings = {cell_idx, 'c'};
        ds1.plot_cell_map(color_mappings, 'enable_class_colors');
        title(sprintf('ds1: Current cell (i=%d) shown in cyan', cell_idx));
        fprintf('  Showing cell map (press any key to return)\n');
        pause;
        datacursormode off;
    end

    % Convenience functions
    %------------------------------------------------------------
    function set_idx1(new_idx)
        idx1 = new_idx;
        idx2 = 1;
    end

    function increment_i
        idx1 = idx1 + 1;
        if idx1 > num_cells1
            cprintf('blue', '  Already at last ds1 cell!\n');
            idx1 = num_cells1;
        end
        idx2 = 1;
    end

    function decrement_i
        idx1 = idx1 - 1;
        idx1 = max(1, idx1);
        idx2 = 1;
    end

end
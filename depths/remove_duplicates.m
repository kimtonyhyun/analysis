function remove_duplicates(ds1, ds2)
% Given cell maps in two planes, look for duplicates by considering
% spatiotemporal alignment of cells, and then "assign" the cell to one of
% the DaySummary's. Note that the labels are mutated in place.

app_name = 'Remove duplicates';

num_cells1 = ds1.num_classified_cells;
num_cells2 = ds2.num_classified_cells;

% Compute all temporal correlations. Results are sorted by the correlation
% value (descending).
corrlist = compute_corrlist(ds1, ds2);
corrlist = sortrows(corrlist, 1, 'ascend'); % Sort by ds1 cell index

% Next, precompute distances between cells using their center-of-mass. Note
% that this computation is performed for all sources, not just those
% classified to be a cell. TODO: Consider allowing for transform?
num_all_cells1 = ds1.num_cells;
num_all_cells2 = ds2.num_cells;

coms1 = cell2mat({ds1.cells.com}); % coms1(:,k) is the COM of the k-th cell
coms2 = cell2mat({ds2.cells.com});

D = zeros(num_all_cells1, num_all_cells2);
for i = 1:num_all_cells1
    for j = 1:num_all_cells2
        D(i,j) = norm(coms1(:,i)-coms2(:,j));
    end
end

hf = figure;

i = 1; % Loops over ds1 cells
j = 1; % Loops over ds2 cells
while (1)
    % Block of corrlist for the i-th cell in ds1
    rows = (1+(i-1)*num_cells2):(i*num_cells2);
    corrlist_i = corrlist(rows,:);
    
    ds1_cell_idx = corrlist_i(j,1);
    ds2_cell_idx = corrlist_i(j,2);
    corr_val = corrlist_i(j,3);
    
    show_corr(ds1, ds1_cell_idx, ds2, ds2_cell_idx, corr_val,...
        'zsc', 'overlay', [], 'zoom_target', 1);
    
    prompt = sprintf('%s (ds1 idx1=%d of %d; ds2 idx2=%d of %d) >> ',...
                      app_name,...
                      i, num_cells1,...
                      j, num_cells2);
    resp = strtrim(input(prompt, 's'));
    
    resp = lower(resp);
    if isempty(resp)
            % Increment ds2 cell index
            j = j + 1;
            j = min(j, num_cells2);
    else
        val = str2double(resp);
        if (~isnan(val)) % Is a number
            if (val == 0)
                % Find the nearest cell
                [~, ds2_cell_idx] = min(D(ds1_cell_idx,:));
                j = find(corrlist_i(:,2)==ds2_cell_idx,1);
            elseif (1 <= val) && (val <= num_cells2)
                j = val;
            end
        else
            switch resp(1)                  
                case {'a', 'c'} % "Assign"
                    val = str2double(resp(2:end));
                    switch val
                        case 1 % Assign to ds1
                            ds2.cells(ds2_cell_idx).label = 'not a cell';
                            fprintf('  Assigned to ds1!\n');
                            
                            i = i + 1;
                            i = min(num_cells1, i);
                            j = 1;
                            
                        case 2 % Assign to ds2
                            ds1.cells(ds1_cell_idx).label = 'not a cell';
                            fprintf('  Assigned to ds2!\n');
                            
                            i = i + 1;
                            i = min(num_cells1, i);
                            j = 1;
                            
                        otherwise
                            cprintf('red', '  Error: Valid inputs are "a1" or "a2"\n');
                    end

                case 'm' % Cell map
                    display_map(ds1_cell_idx);

                case 'n' % Next
                    i = i + 1;
                    i = min(num_cells1, i);
                    j = 1;
                    
                case 'p' % Previous
                    i = i - 1;
                    i = max(1, i);
                    j = 1;

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

end
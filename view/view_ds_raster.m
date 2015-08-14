function view_ds_raster(ds)

% By default, show rasters of classified cells
cell_indices = find(ds.is_cell);
num_cells = length(cell_indices);

% Display settings
cells_per_page = [2 3];
num_cells_per_page = prod(cells_per_page);
num_pages = ceil(num_cells / num_cells_per_page);

page_idx = 1;
while (1)
    clf;
    display_page(page_idx);

    % Ask user for command
    prompt = sprintf('Raster viewer (page %d of %d) >> ', page_idx, num_pages);
    resp = strtrim(input(prompt, 's'));
    
    val = str2double(resp);
    if (~isnan(val)) % Is a number (cell index)
        if ismember(val, cell_indices)
            idx = find(cell_indices==val,1);
            page_idx2 = ceil(idx / num_cells_per_page);
            if (page_idx2 == page_idx)
                
            else % Else, jump to page that contains the cell
                page_idx = page_idx2;
            end
        else
            fprintf('  Sorry, %d is not a valid cell index\n', val);
        end
    else
        resp = lower(resp);
        switch (resp)
            case {'f', 'n'} % Forward/next (page)
                if (page_idx < num_pages)
                    page_idx = page_idx + 1;
                else
                    fprintf('  Already at final page!\n');
                end

            case {'b', 'p'} % Backward/prev (page)
                if (page_idx > 1)
                    page_idx = page_idx - 1;
                else
                    fprintf('  Already at first page!\n');
                end

            case 'q' % Exit
                break;

            otherwise
                fprintf('  Could not parse "%s"\n', resp);

        end % switch
    end
end

    % Helper functions
    %------------------------------------------------------------
    function display_page(page_idx)
        cells_on_page = get_cells_on_page(page_idx);
        for i = 1:length(cells_on_page)
            subplot(cells_per_page(1), cells_per_page(2), i);
        
            cell_idx = cells_on_page(i);
            ds.plot_cell_raster(cell_idx);
            title(sprintf('Cell %d (%s)', cell_idx, ds.cells(cell_idx).label));
        end
        
        function cells_on_page = get_cells_on_page(page_idx)
            index_to_cells = [1+num_cells_per_page*(page_idx-1) num_cells_per_page*page_idx];
            index_to_cells(2) = min(index_to_cells(2), num_cells);
            cells_on_page = cell_indices(index_to_cells(1):index_to_cells(2));
        end % get_cells_on_page
        
    end % display_page

end % view_ds_raster
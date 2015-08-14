function browse_rasters(ds)

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
            fprintf('  Cell %d selected!\n', val);
            
            % Find the page that contains the cell
            page_idx2 = ceil(find(cell_indices==val,1) / num_cells_per_page);
            if (page_idx2 == page_idx)
                % The selected cell is already in the current page, examine
                % the raster in further detail
                display_cell(val);
                fprintf('  Press any key to return to all cells...\n');
                pause;
                
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

    function display_cell(cell_idx)
        % Image of cell
        subplot(3,4,[1 2]);
        imagesc(ds.cells(cell_idx).im);
        axis image;
        title(sprintf('Cell %d (%s)', cell_idx, ds.cells(cell_idx).label));
        
        % Raster of all trials, with correctness
        subplot(3,4,[5 6 9 10]);
        ds.plot_cell_raster(cell_idx);
        title('All trials');
        
        corr_width = 0.025;
        xlim([0 1+corr_width]);
        for i = 1:ds.num_trials
            if ds.trials(i).correct
                corr_color = 'g';
            else
                corr_color = 'r';
            end
            rectangle('Position', [1 i-0.5 corr_width 1],...
                      'FaceColor', corr_color);
        end
        
        % Divide rasters by correctness
        subplot(3,4,3);
        ds.plot_cell_raster(cell_idx, 'correct');
        title('Correct');
        subplot(3,4,4);
        ds.plot_cell_raster(cell_idx, 'incorrect');
        title('Incorrect');
        
        % Divide rasters by start location
        subplot(3,4,7);
        ds.plot_cell_raster(cell_idx, 'start', 'west');
        title('West start');
        subplot(3,4,8);
        ds.plot_cell_raster(cell_idx, 'start', 'east');
        title('East start');
        
        % Divide rasters by end location
        subplot(3,4,11);
        ds.plot_cell_raster(cell_idx, 'end', 'south');
        title('South end');
        subplot(3,4,12);
        ds.plot_cell_raster(cell_idx, 'end', 'north');
        title('North end');
    end % display_cell

end % view_ds_raster
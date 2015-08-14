function view_ds_raster(ds)

% Display settings
cells_per_page = [2 3];
num_cells_per_page = prod(cells_per_page);

cell_indices = find(ds.is_cell);
num_cells = length(cell_indices);
num_pages = ceil(num_cells / num_cells_per_page);

for p = 1:num_pages
    
    cells_on_page = [1+num_cells_per_page*(p-1) num_cells_per_page*p];
    cells_on_page(2) = min(cells_on_page(2), ds.num_cells);
    cells_on_page = cells_on_page(1):cells_on_page(2);
    clf;
    for i = 1:length(cells_on_page)
        subplot(cells_per_page(1), cells_per_page(2), i);
        
        cell_idx = cells_on_page(i);
        ds.plot_cell_raster(cell_idx);
        title(sprintf('Cell %d', cell_idx));
    end
    pause;
    
    
end

    function cells_on_page = get_cells_on_page(page_idx)
        
    end % get_cells_on_page

end % view_ds_raster
function plot_event_map(obj)

processed_cells = [];
processed_counts = [];
for k = find(obj.is_cell)
    ek = obj.cells(k).events;
    if ~isempty(ek)
        if ~strcmp(ek.info.method, 'rejected')
            processed_cells = [processed_cells k]; %#ok<*AGROW>
            processed_counts = [processed_counts size(ek.data, 1)];
        end
    end
end

obj.plot_cell_map_redblue(processed_cells, processed_counts);
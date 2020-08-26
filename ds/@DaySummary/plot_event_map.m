function plot_event_map(obj)

[event_counts, cell_inds] = obj.get_event_counts;
obj.plot_cell_map_redblue(cell_inds, event_counts);
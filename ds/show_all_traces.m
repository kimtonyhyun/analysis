function show_all_traces(ds)
% Show the normalized activity of all classified cells

colors = 'kbr';
num_colors = length(colors);

ind = 1;
for k = 1:ds.num_cells
    if ds.is_cell(k)
        tr = ds.get_trace(k);
        tr_min = min(tr);
        tr_max = max(tr);
        tr = (tr-tr_min)/(tr_max-tr_min);
        
        color = colors(mod(ind,num_colors)+1);
        plot(tr+ind, color);
        hold on;
        
        ind = ind + 1;
    end
end

xlim([1 length(tr)]);
xlabel('Frames');
ylim([0 ind]);
ylabel('Cell index');

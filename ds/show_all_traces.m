function show_all_traces(ds, varargin)
% Show the normalized activity of all classified cells

amp = 1; % Amplitude of each trace (relative to spacing)
for k = 1:length(varargin)
    vararg = varargin{k};
    if ischar(vararg)
        switch lower(vararg)
            case {'amp', 'amplitude'}
                amp = varargin{k+1};
        end
    end
end


colors = 'kbr';
num_colors = length(colors);

ind = 1;
for k = 1:ds.num_cells
    if ds.is_cell(k)
        tr = ds.get_trace(k);
        tr_min = min(tr);
        tr_max = max(tr);
        tr = amp*(tr-tr_min)/(tr_max-tr_min);
        
        color = colors(mod(ind,num_colors)+1);
        plot(tr+ind, color);
        hold on;
        
        ind = ind + 1;
    end
end

xlim([1 length(tr)]);
xlabel('Frames');
ylim([0 (ind-1)+amp]);
ylabel('Cell index');

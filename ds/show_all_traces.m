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

% If no cells are classified, then assume all sources are cells
show_all_cells = (ds.num_classified_cells==0);

colors = 'kbr';
num_colors = length(colors);

ind = 1;
for k = 1:ds.num_cells
    if show_all_cells || ds.is_cell(k)
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

num_frames = length(tr);
y_range = [0.5 (ind-1)+amp+0.5];

if (ds.num_trials > 1) % DaySummary contains trial structure
    trial_start_frames = ds.trial_indices(:,1);
    for k = 1:ds.num_trials
        plot(trial_start_frames(k)*[1 1], y_range, 'k:');
    end
    
    x_ticks = trial_start_frames(1:5:end);
    xticks(x_ticks);
    
    x_tick_labels = 1+5*(0:length(x_ticks)-1);
    xticklabels(x_tick_labels);
    xlabel(sprintf('Trials (%d total)', ds.num_trials));
else
    xlabel(sprintf('Frames (%d total)', num_frames));
end
hold off;

xlim([1 num_frames]);
ylim(y_range);
ylabel(sprintf('Cell index (%d total)', ind-1));
set(gca, 'TickLength', [0 0]);

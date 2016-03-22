function view_trace_in_space(ds, cell_idx)
% Displays the trace of 'cell_idx' over the animal's trajectory in the
% behavior arena

[m, M] = get_trace_min_max(ds, cell_idx);

clf;

subplot(2,4,1);
% HACK: We are getting the common scale for the subsequent raster plots by
% drawing a full raster, which will be "painted over".
ds.plot_cell_raster(cell_idx, 'draw_correct');
raster_scale = get(gca, 'CLim');

ds.plot_cell_raster(cell_idx, 'start', 'west', 'draw_correct');
set(gca, 'CLim', raster_scale);
title('West start');

subplot(2,4,8);
ds.plot_cell_raster(cell_idx, 'start', 'east', 'draw_correct');
set(gca, 'CLim', raster_scale);
title('East start');

subplot(2,4,5);
ds.plot_cell_raster(cell_idx, 'end', 'south', 'draw_correct');
set(gca, 'CLim', raster_scale);
title('South end');

subplot(2,4,4);
ds.plot_cell_raster(cell_idx, 'end', 'north', 'draw_correct');
set(gca, 'CLim', raster_scale);
title('North end'); 

subplot(2,4,[2 3 6 7]);
bg_image = ds.behavior_ref_img;
bg_image = cat(3, bg_image, bg_image, bg_image);
imagesc(bg_image);
axis image;
hold on;
title(sprintf('Cell %d (%s)', cell_idx, ds.cells(cell_idx).label));

for trial_idx = 1:ds.num_trials
    centroids = ds.trials(trial_idx).centroids;
    trace = ds.get_trace(cell_idx, trial_idx);
    trace = 255*(trace-m)/(M-m); % Scale trace to image colormap
    cline(centroids(:,1), centroids(:,2), [], trace);
end
hold off;

end % view_trace_in_space

function [global_min, global_max] = get_trace_min_max(ds, cell_idx)
    % Read all trials from DaySummary, so that we can apply a common 
    % scaling to the trace
    global_min = Inf; % Global min
    global_max = -Inf; % Global max
    
    for trial_idx = 1:ds.num_trials
        trace = ds.get_trace(cell_idx, trial_idx);
        m = min(trace);
        if m < global_min
            global_min = m;
        end
        M = max(trace);
        if M > global_max
            global_max = M;
        end
    end
end
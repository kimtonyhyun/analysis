function draw_path_subrasters(ds, cell_idx, raster_scale)
% Used by 'view_raster_touch'

sd = ds.switchdata;

% First column is for CONSTANT path
subplot(2,4,3);
[~, trial_inds] = ds.plot_cell_raster(cell_idx,...
    'range', sd.pre_switch_trials,...
    'start', sd.constant_path_start,...
    'correct', 'draw_correct');
set(gca, 'CLim', raster_scale);
title(sprintf('Constant path (%s-start)', sd.constant_path_start));
n_trials = length(trial_inds);
ylabel(sprintf('PRE-switch (%d)', n_trials));
yticks(1:n_trials);
yticklabels(num2cell(trial_inds));

subplot(2,4,7);
[~, trial_inds] = ds.plot_cell_raster(cell_idx,...
    'range', sd.post_switch_trials,...
    'start', sd.constant_path_start,...
    'correct', 'draw_correct');
set(gca, 'CLim', raster_scale);
title(sprintf('Constant path (%s-start)', sd.constant_path_start));
n_trials = length(trial_inds);
ylabel(sprintf('POST-switch (%d)', n_trials));
yticks(1:n_trials);
yticklabels(num2cell(trial_inds));

% Second column is for CHANGING path
subplot(2,4,4);
[~, trial_inds] = ds.plot_cell_raster(cell_idx,...
    'range', sd.pre_switch_trials,...
    'start', sd.changing_path_start,...
    'correct', 'draw_correct');
set(gca, 'CLim', raster_scale);
title(sprintf('Changing path (%s-start)', sd.changing_path_start));
n_trials = length(trial_inds);
ylabel(sprintf('PRE-switch (%d)', n_trials));
yticks(1:n_trials);
yticklabels(num2cell(trial_inds));

subplot(2,4,8);
[~, trial_inds] = ds.plot_cell_raster(cell_idx,...
    'range', sd.post_switch_trials,...
    'start', sd.changing_path_start,...
    'correct', 'draw_correct');
set(gca, 'CLim', raster_scale);
title(sprintf('Changing path (%s-start)', sd.changing_path_start));
n_trials = length(trial_inds);
ylabel(sprintf('POST-switch (%d)', n_trials));
yticks(1:n_trials);
yticklabels(num2cell(trial_inds));

end % draw_path_subrasters
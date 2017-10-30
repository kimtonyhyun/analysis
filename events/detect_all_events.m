function events = detect_all_events(ds)

if ds.full_num_frames ~= ds.trial_indices(end,end)
    fprintf('  Warning: Event detection should be performed on all trials, including probes!\n');
end

events = cell(ds.num_cells, 1);
for k = 1:ds.num_cells
    if (k==1) || (mod(k,50)==0)
        fprintf('%s: At cell %d of %d...\n', datestr(now), k, ds.num_cells);
    end
    events{k} = detect_events(ds, k, 'noprompt');
end

events = cell2mat(events); % Convert cell to array of structs

% Save to file
timestamp = datestr(now, 'yymmdd-HHMMSS');
event_savename = sprintf('events_%s.mat', timestamp);
save(event_savename, 'events', '-v7.3');
function events = detect_all_events(ds)

events = cell(ds.num_cells, 1);
for k = 1:ds.num_cells
    if (k==1) || (mod(k,50)==0)
        fprintf('%s: At cell %d of %d...\n', datestr(now), k, ds.num_cells);
    end
    events{k} = detect_events(ds, k, 'noprompt');
end

events = cell2mat(events); % Convert cell to array of structs
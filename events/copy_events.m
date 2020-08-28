function copy_events(ds_source, ds_target, match)
% Copy event parameters from source DaySummary to target DaySummary. If no  
% external 'match' provided, copy_labels will compute the matching
% but assuming no translation needed.
%

if (nargin == 2)
    match = run_alignment(ds_source, ds_target, 'notrans', 'noprompt');
end

num_events_copied = 0;

ds_target.reset_labels();
for k = 1:ds_source.num_cells
    m = match{k};
    if ~isempty(m)
        source_events = ds_source.cells(k).events;
        if ~isempty(source_events)
            k2 = m(1,1);
            ds_target.cells(k2).events = ds_source.cells(k).events;
            num_events_copied = num_events_copied + 1;
        end
    end
end
fprintf('  %s: Copied %d events from "%s" to "%s"\n',...
    datestr(now), num_events_copied, inputname(1), inputname(2));
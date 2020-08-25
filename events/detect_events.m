function detect_events(ds, fps, varargin)
% A convenience wrapper around 'detect_events_interactively.m'. Mirrors the
% design of 'classify_cells.m'

M = [];
for j = 1:length(varargin)
    vararg = varargin{j};
    if ischar(vararg)
        switch lower(vararg)
            case {'m', 'movie'}
                M = varargin{j+1};
        end
    end
end

hfig = figure;
if ~isempty(M)
    movie_clim = compute_movie_scale(M);
end

cell_idx = 1;
while (1)
    if isempty(ds.cells(cell_idx).events)
        event_str = 'unprocessed';
    else
        num_events = size(ds.cells(cell_idx).events.data, 1);
        if num_events == 1
            event_str = '1 event';
        else
            event_str = sprintf('%d events', num_events);
        end
    end
    prompt = sprintf('Detector (%d/%d, %s, %s) >> ',...
        cell_idx, ds.num_cells, ds.cells(cell_idx).label, event_str);
    
    resp = strtrim(input(prompt, 's'));
    val = str2double(resp);

    if (~isnan(val) && isreal(val)) % Is a number
        if ((1 <= val) && (val <= ds.num_cells))
            cell_idx = val;
        else
            fprintf('  Sorry, %d is not a valid cell index\n', val);
        end
    else % Not a number
        resp = lower(resp);
        if isempty(resp)
            resp = 'n';
        end
        
        switch (resp(1))
            case 'q' % "quit"
                break;
                
            case {'d', 'c'} % Detect events
                detect_events_interactively(ds, cell_idx, fps, 'hfig', hfig);
                
            case 'm' % Detect events with movie
                if ~isempty(M)
                    detect_events_interactively(ds, cell_idx, fps, 'hfig', hfig,...
                        'movie', M, 'movie_clim', movie_clim);
                else
                    fprintf('  Movie not loaded!\n');
                end
                
            case 'n' % Next cell
                go_to_next_cell();
                
            otherwise
                fprintf('  Sorry, could not parse "%s"\n', resp);
        end
    end
end % Main interaction loop

    function go_to_next_cell()
        no_events = cellfun(@isempty, {ds.cells.events});
        next_idx = find_next_cell_to_process(cell_idx, ds.is_cell & no_events);
        
        if isempty(next_idx)
            cprintf([0 0.5 0], '  All cells have been classified!\n');
        else
            cell_idx = next_idx;
        end
    end

end % detect_events
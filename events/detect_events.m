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

% At the end of 'detect_events', will save the new event data to file if it
% has changed from its initial values.
initial_ds_events = {ds.cells.events};

events_rejected = struct('info', [], 'data', []);
events_rejected.info = struct('method', 'rejected');

cell_idx = 1;
while (1)
    display_trace();
    
    e = ds.cells(cell_idx).events;
    if isempty(e)
        event_str = 'unprocessed';
    else
        switch e.info.method
            case 'rejected'
                event_str = 'rejected';
            otherwise
                num_events = size(e.data, 1);
                if num_events == 1
                    event_str = '1 event';
                else
                    event_str = sprintf('%d events', num_events);
                end
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
            go_to_next_cell()
        else
            switch (resp(1))
                case 'q' % "quit"
                    if ~isequal(initial_ds_events, {ds.cells.events})
                        savename = ds.save_events;
                        fprintf('  Saved events data to "%s"\n', savename);
                    else
                        fprintf('  No changes to events data\n');
                    end
                    close(hfig);
                    break;

                case 's' % save
                    savename = ds.save_events;
                    fprintf('  Saved events data to "%s"\n', savename);
                    
                case 'l' % Load previous classification
                    [file, path] = uigetfile('events_*.mat', 'Select existing events file');
                    if (file)
                        full_file = fullfile(path, file);
                        ds.load_events(full_file);
                    end
                    
                case {'d', 'c'} % Detect events
                    events = detect_events_interactively(ds, cell_idx, fps, 'hfig', hfig);
                    
                    % Write latest events to DaySummary?
                    ds_events = ds.cells(cell_idx).events;
                    if ~isequal(ds_events, events)
                        if isempty(ds_events)
                            fprintf('  Save event results to DaySummary? (Y/n) ');
                            default = 'y';
                        else
                            cprintf('red', '  Overwrite existing event results in DaySummary? (y/N) ');
                            default = 'n';
                        end
                        resp = strtrim(input('>> ', 's'));
                        if isempty(resp)
                            resp = default;
                        end
                        switch lower(resp)
                            case 'y'
                                ds.cells(cell_idx).events = events;
                        end
                    end

                case 'm' % Detect events with movie
                    if ~isempty(M)
                        detect_events_interactively(ds, cell_idx, fps, 'hfig', hfig,...
                            'movie', M, 'movie_clim', movie_clim);
                    else
                        fprintf('  Movie not provided!\n');
                    end

                case {'n', 'r'} % Reject cell for event detection (e.g. if too noisy)
                    ds.cells(cell_idx).events = events_rejected;
                    fprintf('  Cell %d rejected for event detection\n', cell_idx);
                    go_to_next_cell();
                    
                case 'p'
                    if cell_idx > 1
                        cell_idx = cell_idx - 1;
                    end

                otherwise
                    fprintf('  Sorry, could not parse "%s"\n', resp);
            end
        end
    end
end % Main interaction loop

    function display_trace()
        clf(hfig, 'reset');
        subplot(311);
        ds.plot_trace(cell_idx);
        title(sprintf('Cell %d of %d', cell_idx, ds.num_cells));
        
        subplot(3,2,[3 5]);
        ds.plot_event_map;
        hold on;
        to_be_processed = find(ds.is_cell & cellfun(@isempty, {ds.cells.events}));
        for ind = to_be_processed
            bd = ds.cells(ind).boundary;
            plot3(bd(:,1), bd(:,2),...
                  ind*ones(size(bd,1),1),... % Allows display of cell number by clicking boundary (as in plot_cell_map)
                  'w:');
        end
        com = ds.cells(cell_idx).com;
        plot(com(1), com(2), 'c.', 'MarkerSize', 24);
        hold off;
    end % display_trace

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
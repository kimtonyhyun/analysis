function classify_cells(ds, M, varargin)
% Perform manual classification of candidate filter/trace pairs

% Note: Default options configured for the PlusMaze experiment!
show_raster = true;
fps = 10;

for i = 1:length(varargin)
    vararg = varargin{i};
    if ischar(vararg)
        switch lower(vararg)
            case 'fps'
                fps = varargin{i+1};
            case 'noraster' % Use with non-PlusMaze datasets
                show_raster = false;
        end
    end
end

% Initial processing of movie
max_proj = max(M,[],3);
movie_clim = compute_movie_scale(M);
fprintf('  %s: Movie will be displayed with fixed CLim = [%.3f %.3f]...\n',...
    datestr(now), movie_clim(1), movie_clim(2));
fprintf('  %s: FPS is %.1f...\n', datestr(now), fps);

% Load filter/trace pairs to be classified
num_candidates = ds.num_cells;

assert(size(M,3) == ds.full_num_frames,...
       'Number of frames in movie does not match that in DaySummary!');

% Begin classification
%------------------------------------------------------------
output_name = sprintf('class_%s.txt', datestr(now, 'yymmdd-HHMMSS'));

hfig = figure;

cell_idx = 1;
prev_cell_idx = 1;

while (cell_idx <= num_candidates)
    if show_raster
        display_candidate_rasters(cell_idx);
    else
        display_candidate(cell_idx);
    end
    
    % Ask the user to classify the cell candidate
    prompt = sprintf('Classifier (%d/%d, "%s") >> ', ...
                        cell_idx, num_candidates, ds.cells(cell_idx).label);
    resp = strtrim(input(prompt, 's'));
    
    val = str2double(resp);
    if (~isnan(val)) % Is a number. Check if it is a valid index and jump to it
        if ((1 <= val) && (val <= num_candidates))
            prev_cell_idx = cell_idx;
            cell_idx = val;
        else
            fprintf('  Sorry, %d is not a valid cell index\n', val);
        end
    else

        % Apply syntactic sugar
        switch (resp)
            case 'C'
                resp = 'c!';
            case 'N'
                resp = 'n';
        end

        resp = lower(resp);
        switch (resp)
            % Classication options
            %------------------------------------------------------------
            case {'p', 'c'} % Cell
                [~, movie_clim] = view_cell_interactively(ds, cell_idx,...
                                    M, fps, movie_clim);
                resp2 = input(sprintf('  Confirm classification ("%s") >> ', resp), 's');
                resp2 = lower(strtrim(resp2));
                if (strcmp(resp, resp2)) % Confirmed
                    set_label(cell_idx, resp2);
                    go_to_next_unlabeled_cell();
                end

            case {'p!', 'c!', 'n'} % Classify without viewing trace
                set_label(cell_idx, resp(1));
                go_to_next_unlabeled_cell();
                
            % Application options
            %------------------------------------------------------------
            case ''  % Go to next unlabeled cell candidate, loop at end
                go_to_next_unlabeled_cell();
            case '-' % Jump to previously viewed cell
                temp_idx = cell_idx;
                cell_idx = prev_cell_idx;
                prev_cell_idx = temp_idx;
            case 'm' % View cell map
                display_map();
            case 'q' % Exit
                close(hfig);
                break;
            case 's' % Save classification
                ds.save_class(output_name);
                fprintf('  Saved classification result to %s\n', output_name);
            case 'l' % Load previous classification
                [file, path] = uigetfile('*.txt', 'Select existing classification');
                if (file)
                    full_file = fullfile(path, file);
                    ds.load_class(full_file);
                end
            case 't' % "Take" screenshot
                screenshot_name = sprintf('cell%03d.png', cell_idx);
                print('-dpng', screenshot_name);
                fprintf('  Plot saved to %s\n', screenshot_name);
            otherwise
                fprintf('  Could not parse "%s"\n', resp);
        end
    end
end

% Save at end!
ds.save_class(output_name);

    % Auxiliary functions
    %------------------------------------------------------------
    function display_candidate_rasters(cell_idx)
        clf;
        subplot(3,2,[1 2]);
        ds.plot_trace(cell_idx);
        title(sprintf('Candidate %d of %d', cell_idx, num_candidates));
        
        subplot(3,2,[3 5]);
        ds.plot_superposed_trials(cell_idx);

        subplot(3,2,[4 6]);
        ds.plot_cell_raster(cell_idx, 'draw_correct');
    end % display_candidate_rasters

    function display_candidate(cell_idx)
        clf;
        subplot(3,1,1);
        ds.plot_trace(cell_idx);
        title(sprintf('Candidate %d of %d', cell_idx, num_candidates));
        
        % Plot cell filter on top of max projection image
        subplot(3,1,[2 3]);
        imagesc(max_proj);
        axis image;
        colormap gray;
        hold on;
        
        filter_threshold = 0.3;
        filter = ds.cells(cell_idx).im;
        boundaries = compute_ic_boundary(filter, filter_threshold);
        for j = 1:length(boundaries)
            boundary = boundaries{j};
            plot(boundary(:,1), boundary(:,2), 'c', 'LineWidth', 2);
        end
        
        COM = ds.cells(cell_idx).com;
        plot(COM(1), COM(2), 'b.');
        
        % Draw nearest neighbors
        num_neighbors_to_draw = min(20, ds.num_cells-1);
        other_cells = ds.get_nearest_sources(cell_idx, num_neighbors_to_draw);
        for oc_idx = other_cells
            oc = ds.cells(oc_idx);
            if isempty(oc.label)
                color = 'w';
            else
                if ds.is_cell(oc_idx)
                    color = 'g';
                else
                    color = 'r';
                end
            end
            plot(oc.boundary(:,1), oc.boundary(:,2), color);
        end
        hold off;
        
        [height, width, ~] = size(M);
        zoom_half_width = min([width, height])/10;
        x_range = COM(1)+zoom_half_width*[-1 1];
        y_range = COM(2)+zoom_half_width*[-1 1];
        xlim(x_range);
        ylim(y_range);
    end

    function display_map()
        clf;
        color_mappings = {cell_idx, 'c'};
        ds.plot_cell_map(color_mappings, 'enable_class_colors');
        title(sprintf('Current cell (ID=%d) shown in cyan', cell_idx));
        fprintf('  Showing cell map (press any key to return)\n');
        pause;
        datacursormode off;
    end % display_map
    
    function go_to_next_unlabeled_cell()
        labels = ds.get_class;
        unlabeled = cellfun(@isempty, labels);
        unlabeled = circshift(unlabeled, -cell_idx);
        search_offset = find(unlabeled, 1);
        if isempty(search_offset)
            fprintf('  All cells have been classified!\n');
        else
            prev_cell_idx = cell_idx;
            cell_idx = mod(cell_idx+search_offset-1, num_candidates) + 1;
        end
    end
        
    function set_label(cell_idx, label)
        switch label
            case 'p'
                full_label = 'phase-sensitive cell';
            case 'c'
                full_label = 'cell';
            case 'n'
                full_label = 'not a cell';
        end
        ds.cells(cell_idx).label = full_label;
        fprintf('  Candidate %d classified as %s\n', cell_idx, full_label);
    end % set_label
end % classify_cells

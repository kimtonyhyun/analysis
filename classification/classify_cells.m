function classify_cells(ds, M)
% Perform manual classification of candidate filter/trace pairs

% Compute a common scaling for the movie
movie_clim = compute_movie_scale(M);
fprintf('  %s: Movie will be displayed with fixed CLim = [%.3f %.3f]...\n',...
    datestr(now), movie_clim(1), movie_clim(2));

fps = 10; % FIXME

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
    display_candidate(cell_idx);
    
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
    function display_candidate(cell_idx)
        clf;
        subplot(3,4,[1 2 3]);
        ds.plot_trace(cell_idx);
        title(sprintf('Candidate %d of %d', cell_idx, num_candidates));

        subplot(3,4,4);
        display_neighborhood(cell_idx);
        
        subplot(3,2,[3 5]);
        ds.plot_superposed_trials(cell_idx);

        subplot(3,2,[4 6]);
        ds.plot_cell_raster(cell_idx, 'draw_correct');
        
        function display_neighborhood(cell_idx)
            % Display the current cell
            cell_img = ds.cells(cell_idx).im;
            [height, width] = size(cell_img);
            zoom_half_width = min([height, width])/15;
            
            imagesc(cell_img);
            axis image;
            xlabel('x [px]');
            ylabel('y [px]');
            com = ds.cells(cell_idx).com;
            x_range = com(1) + zoom_half_width*[-1 1];
            x_range(1) = max(1, x_range(1)); x_range(2) = min(width, x_range(2));
            y_range = com(2) + zoom_half_width*[-1 1];
            y_range(1) = max(1, y_range(1)); y_range(2) = min(height, y_range(2));
            xlim(x_range);
            ylim(y_range);
            hold on;
            
            % Display neighbors
            num_neighbors_to_show = 5;
            neighbor_indices = ds.get_nearest_sources(cell_idx, num_neighbors_to_show);
            for k = 1:num_neighbors_to_show
                neighbor_idx = neighbor_indices(k);
                boundary = ds.cells(neighbor_idx).boundary;
                com = ds.cells(neighbor_idx).com;
                
                color = 'w';
                plot(boundary(:,1), boundary(:,2), color);
                text(com(1), com(2), num2str(neighbor_idx),...
                     'HorizontalAlignment', 'center',...
                     'Color', color,...
                     'Clipping', 'on');
            end
        end % Display neighborhood
    end % display_candidate

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

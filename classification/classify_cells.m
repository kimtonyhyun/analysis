function poi = classify_cells(ds, M, varargin)
% Perform manual classification of candidate filter/trace pairs

fps = 10; % Not super critical

% Default parameters for "view_cell_interactively"
state.show_map = true;
state.show_neighbors = false;
state.threshold_scale = 0.5;
state.points_of_interest = [];

% Will be set below
state.movie_clim = []; % Used within view_cell_interactively
state.ref_image_clim = [];
state.fig_handle = [];

% Custom "subplot" command that leaves less unusued space between panels
sp = @(m,n,p) subtightplot(m, n, p, 0.05, 0.05, 0.05); % Gap, Margin-X, Margin-Y

for i = 1:length(varargin)
    vararg = varargin{i};
    if ischar(vararg)
        switch lower(vararg)
            case 'fps'
                fps = varargin{i+1};
            case 'poi' % Existing points of interest
                state.points_of_interest = varargin{i+1};
        end
    end
end

num_frames = size(M,3);

% Initial processing of movie
if num_frames > 1
    if ~(num_frames == ds.full_num_frames)
        cprintf('blue',...
           'Number of frames in movie (%d) does not match trace length (%d) in DaySummary!\n',...
           size(M,3), ds.full_num_frames);
    end
    
    ref_image = max(M,[],3);
    ref_image_title = 'Movie max projection';
    
    state.movie_clim = compute_movie_scale(M);
    fprintf('  %s: Movie will be displayed with fixed CLim = [%.3f %.3f]...\n',...
        datestr(now), state.movie_clim(1), state.movie_clim(2));    
    fprintf('  %s: FPS is %.1f...\n', datestr(now), fps);
else
    % Single image
    cprintf('blue', 'Classifying filters with respect to an image!\n');
    ref_image = M; % The provided image itself
    ref_image_title = 'Reference image';
end
state.ref_image_clim = compute_movie_scale(ref_image);

% Load filter/trace pairs to be classified
num_candidates = ds.num_cells;

% Begin classification
%------------------------------------------------------------
timestamp = datestr(now, 'yymmdd-HHMMSS');
output_name = sprintf('class_%s.txt', timestamp);

state.fig_handle = figure;

cell_idx = 1;
prev_cell_idx = 1;

while (cell_idx <= num_candidates)
    if num_frames > 1
        % Classify cells with respect to calcium movie
        if ds.num_trials > 1
            display_candidate_rasters(cell_idx);
        else
            display_candidate(cell_idx, false);
        end
    else
        % Classify cells with respect to a provided image
        display_candidate(cell_idx, true);
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

        % Some syntactic sugar
        switch (resp)
            case 'C'
                resp = 'c!';
            case 'N'
                resp = 'n';
            case 'Q'
                resp = 'q!';
        end

        resp = lower(resp);
        switch (resp)
            % Classication options
            %------------------------------------------------------------
            case 'c' % Cell
                if num_frames > 1 % Actual movie
                    [resp2, state] = view_cell_interactively(ds, cell_idx, M, fps, state);
                    switch resp2
                        % For these "exit codes", label the candidate and move on
                        case {'c', '+'}
                            set_label(cell_idx, 'c');
                            go_to_next_unlabeled_cell();

                        case {'n', '-'}
                            set_label(cell_idx, 'n');
                            go_to_next_unlabeled_cell();

                        % "Classic" behavior
                        otherwise
                            resp2 = input(sprintf('  Confirm classification ("%s") >> ', resp), 's');
                            resp2 = lower(strtrim(resp2));
                            if (strcmp(resp, resp2)) % Confirmed
                                set_label(cell_idx, resp2);
                                go_to_next_unlabeled_cell();
                            end
                    end
                else % Single image
                    set_label(cell_idx, 'c');
                    go_to_next_unlabeled_cell();
                end

            case 'c!' % Classify without viewing trace
                set_label(cell_idx, 'c');
                go_to_next_unlabeled_cell();
                
            case 'n' % Classify without viewing trace
                set_label(cell_idx, 'n');
                go_to_next_unlabeled_cell();
                
            % Application options
            %------------------------------------------------------------
            case 'h' % Higher contrast
                c_range = diff(state.ref_image_clim);
                state.ref_image_clim = state.ref_image_clim + c_range*[0.1 -0.1];
            case 'l' % Lower contrast
                c_range = diff(state.ref_image_clim);
                state.ref_image_clim = state.ref_image_clim + c_range*[-0.1 0.1];

            case ''  % Go to next unlabeled cell candidate, loop at end
                go_to_next_unlabeled_cell();
            case {'p', '-'} % Jump to previously viewed cell
                temp_idx = cell_idx;
                cell_idx = prev_cell_idx;
                prev_cell_idx = temp_idx;
            case 'm' % View cell map
                display_map();
            case 'q' % Exit
                % Save results
                ds.save_class(output_name);
                poi = state.points_of_interest;
                if ~isempty(poi)
                    num_points = size(poi,1);
                    poi_savename = sprintf('poi_%s.mat', timestamp);
                    save(poi_savename, 'poi');
                    if (num_points == 1)
                        pt_str = 'point';
                    else
                        pt_str = 'points';
                    end
                    fprintf('  Saved %d %s of interest to "%s"\n',...
                        num_points, pt_str, poi_savename);
                end
                
                close(state.fig_handle);
                break;
            case 'q!' % Exit without saving
                close(state.fig_handle);
                break;
            case 's' % Save classification
                ds.save_class(output_name);
                fprintf('  Saved classification result to %s\n', output_name);
            case 'l' % Load previous classification
                [file, path] = uigetfile('class_*.txt', 'Select existing classification');
                if (file)
                    full_file = fullfile(path, file);
                    ds.load_class(full_file);
                end
            case 't' % "Take" screenshot
                screenshot_name = sprintf('cell%03d.png', cell_idx);
                print('-dpng', screenshot_name);
                fprintf('  Plot saved to %s\n', screenshot_name);
            case 'f'
                
            otherwise
                fprintf('  Could not parse "%s"\n', resp);
                
        end
    end
end

    % Auxiliary functions
    %------------------------------------------------------------
    function display_candidate_rasters(cell_idx)
        clf;
        sp(3,2,[1 2]);
        ds.plot_trace(cell_idx);
        title(sprintf('Source %d of %d', cell_idx, num_candidates));
        
        sp(3,2,[3 5]);
        ds.plot_superposed_trials(cell_idx);

        sp(3,2,[4 6]);
        ds.plot_cell_raster(cell_idx, 'draw_correct');
    end % display_candidate_rasters

    function display_candidate(cell_idx, use_ref_image)
        clf;
        sp(3,1,1);
        ds.plot_trace(cell_idx);
        title(sprintf('Source %d of %d', cell_idx, num_candidates));
        
        COM = ds.cells(cell_idx).com;

        sp(3,2,[3 5]);
        imagesc(ref_image, state.ref_image_clim);
        title(ref_image_title);
        colormap gray;
        axis image;
        hold on;
        
        filter_threshold = 0.3;
        filter = ds.cells(cell_idx).im;
        boundaries = compute_ic_boundary(filter, filter_threshold);
        for j = 1:length(boundaries)
            boundary = boundaries{j};
            plot(boundary(:,1), boundary(:,2), 'c', 'LineWidth', 2, 'HitTest', 'off');
        end
        plot(COM(1), COM(2), 'b.', 'HitTest', 'off');
        
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
            plot(oc.boundary(:,1), oc.boundary(:,2), color, 'HitTest', 'off');
            text(oc.com(1), oc.com(2), num2str(oc_idx),...
                 'HitTest', 'off',...
                 'Clipping', 'on',...
                 'HorizontalAlignment', 'center',...
                 'Color', 'w',...
                 'FontWeight', 'bold');
        end
        
        % Draw points of interest
        pois = state.points_of_interest;
        if ~isempty(pois)
            plot(pois(:,1), pois(:,2), 'y*');
        end
        hold off;
        
        zoom_half_width = 50;
        x_range = COM(1)+zoom_half_width*[-1 1];
        y_range = COM(2)+zoom_half_width*[-1 1];
        xlim(x_range);
        ylim(y_range);
        
        sp(3,2,[4 6]);
        if use_ref_image
            C = imfuse(ref_image, ds.cells(cell_idx).im,...
                    'falsecolor', 'ColorChannels', [1 2 0]);
            image(C);
            title('Spatial filter (green); Ref image (red)');
        else
            imagesc(ds.cells(cell_idx).im);
            title('Spatial filter');
            colormap gray;
        end
        axis image;
        xlim(x_range);
        ylim(y_range);
        set(gca, 'YTick', []);
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
        unlabeled = cellfun(@isempty, ds.get_class);
        next_idx = find_next_cell_to_process(cell_idx, unlabeled);
        
        if isempty(next_idx)
            cprintf([0 0.5 0], '  All cells have been classified!\n');
            prev_cell_idx = cell_idx;
            cell_idx = cell_idx + 1;
            if cell_idx > ds.num_cells
                cell_idx = 1;
            end
        else
            prev_cell_idx = cell_idx;
            cell_idx = next_idx;
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

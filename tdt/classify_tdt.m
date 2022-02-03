function classify_tdt(ds, ref_image)
% Originally based on 'classify_cells'

% Will be set below
state.ref_image_clim = [];
state.fig_handle = [];

% Custom "subplot" command that leaves less unusued space between panels
sp = @(m,n,p) subtightplot(m, n, p, 0.05, 0.05, 0.05); % Gap, Margin-X, Margin-Y

% for i = 1:length(varargin)
%     vararg = varargin{i};
%     if ischar(vararg)
%         switch lower(vararg)
%             case 'fps'
%                 fps = varargin{i+1};
%             case 'poi' % Existing points of interest
%                 state.points_of_interest = varargin{i+1};
%         end
%     end
% end


% Initial processing of movie
ref_image_title = 'Reference image';

state.ref_image_clim = compute_movie_scale(ref_image);

% Load filter/trace pairs to be classified
num_candidates = ds.num_cells;

% Begin classification
%------------------------------------------------------------
timestamp = datestr(now, 'yymmdd-HHMMSS');
output_name = sprintf('tdt_%s.mat', timestamp);

state.fig_handle = figure;

cell_idx = 1;
prev_cell_idx = 1;

while (cell_idx <= num_candidates)

    display_candidate(cell_idx);

    % Ask the user to classify the cell candidate
    prompt = sprintf('tdT classifier (%d/%d, "%s") >> ', ...
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
            case {'c', 'c!'} % Cell
                set_label(cell_idx, 'c');
                go_to_next_unlabeled_cell();

            case 'n' % Classify without viewing trace
                set_label(cell_idx, 'n');
                go_to_next_unlabeled_cell();
                
            % Application options
            %------------------------------------------------------------
            case 'hc' % Higher contrast
                c_range = diff(state.ref_image_clim);
                state.ref_image_clim = state.ref_image_clim + c_range*[0.1 -0.1];
            case 'lc' % Lower contrast
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
                
                close(state.fig_handle);
                break;
            case 'q!' % Exit without saving
                close(state.fig_handle);
                break;

            case 's' % Save classification
                save_tdt(ds, output_name);
            case 'l' % Load previous classification
                [file, path] = uigetfile('tdt_*.mat', 'Select existing classification');
                if (file)
                    full_file = fullfile(path, file);
                    load_tdt(ds, full_file);
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

    % Auxiliary functions
    %------------------------------------------------------------
    function display_candidate(cell_idx)
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
                switch(oc.label)
                    case 'positive'
                        color = [0.85 0.325 0.098]; % Orange
                    case 'negative'
                        color = [0 0.4470 0.7410]; % Blue
                end
            end
            plot(oc.boundary(:,1), oc.boundary(:,2), 'Color', color, 'HitTest', 'off');
            text(oc.com(1), oc.com(2), num2str(oc_idx),...
                 'HitTest', 'off',...
                 'Clipping', 'on',...
                 'HorizontalAlignment', 'center',...
                 'Color', 'w',...
                 'FontWeight', 'bold');
        end
        hold off;
        
        zoom_half_width = 50;
        x_range = COM(1)+zoom_half_width*[-1 1];
        y_range = COM(2)+zoom_half_width*[-1 1];
        xlim(x_range);
        ylim(y_range);
        
        sp(3,2,[4 6]);
        C = imfuse(ref_image, ds.cells(cell_idx).im,...
                   'falsecolor', 'ColorChannels', [1 2 0]);
        image(C);
        title('Spatial filter (green); Ref image (red)');
        axis image;
        xlim(x_range);
        ylim(y_range);
        set(gca, 'YTick', []);
    end

    function display_map()
        clf;
        [pos, neg] = collect_tdt_labels(ds);
        pos = setdiff(pos, cell_idx);
        neg = setdiff(neg, cell_idx);
        unlabeled = setdiff(1:ds.num_cells, [pos neg cell_idx]);

        color_grouping = {cell_idx, 'c';
                          pos, [0.85 0.325 0.098]; % Orange
                          neg, [0 0.447 0.741]; % Blue
                          unlabeled, 'w'};
        ds.plot_cell_map(color_grouping);
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
            case 'c'
                full_label = 'positive';
            case 'n'
                full_label = 'negative';
        end
        ds.cells(cell_idx).label = full_label;
        fprintf('  Candidate %d classified as %s\n', cell_idx, full_label);
    end % set_label
end % classify_cells

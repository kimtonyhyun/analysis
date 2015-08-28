% Summary of PlusMaze data for a single day.
%
% Inputs:
%   plusmaze_txt: Name of the PlusMaze output text file.
%   rec_dir: Directory containing 
%       - Filters and traces in a "rec_*.mat" file (required)
%       - Classification results in a "class_*.txt" file (optional)
%
% Output:
%   DaySummary object
%
% Example usage:
%   ds = DaySummary('c11m1d12_ti2.txt', 'rec001');
%
classdef DaySummary < handle
    properties
        cells
        trials
        
        num_cells
        num_trials
        
        trial_indices
    end
    
    methods
        function obj = DaySummary(plusmaze_txt, rec_dir, varargin)
            % Handle optional input
            exclude_probe_trials = 0;
            for k = 1:length(varargin)
                if ischar(varargin{k})
                    switch lower(varargin{k})
                        case 'excludeprobe'
                            exclude_probe_trials = 1;
                    end
                end
            end
            
            % Load data
            %------------------------------------------------------------
            data_source = get_most_recent_file(rec_dir, 'rec_*.mat');
            data = load(data_source);
            obj.num_cells = data.info.num_pairs;
            fprintf('  %s: Loaded data from %s\n', datestr(now), data_source);
            
            % Parse trial data
            %   TODO: Bring in centroids corresponding to mouse position
            %------------------------------------------------------------
            [trial_indices, loc_info, trial_durations] =...
                parse_plusmaze(plusmaze_txt); %#ok<*PROP>
            
            % Check that the length of traces is consistent with the table
            % of trial indices
            assert(size(data.traces,1) == trial_indices(end,end),...
                'Error: Length of traces does not match trial index table!');
            
            if (exclude_probe_trials)
                is_probe = strcmp(loc_info(:,1), 'north') | ...
                           strcmp(loc_info(:,1), 'south');
                       
                trial_indices = trial_indices(~is_probe,:);
                loc_info = loc_info(~is_probe,:);
                trial_durations = trial_durations(~is_probe);
            end
            
            num_trials = size(trial_indices, 1);
            turns = cell(num_trials, 1);
            traces = cell(num_trials, 1);
            for k = 1:num_trials
                trial_frames = trial_indices(k,1):...
                               trial_indices(k,end);
                traces{k} = data.traces(trial_frames, :)';
                turns{k} = obj.compute_turn(loc_info{k,1}, loc_info{k,3});
            end
            
            obj.num_trials = num_trials;
            obj.trial_indices = trial_indices;
            obj.trials = struct(...
                'start', loc_info(:,1),...
                'goal',  loc_info(:,2),...
                'end',   loc_info(:,3),...
                'correct', cellfun(@strcmp, loc_info(:,2), loc_info(:,3), 'UniformOutput', false),...
                'turn',  turns,...
                'time',  num2cell(trial_durations),...
                'traces', traces);
            
            % Parse cell data
            %------------------------------------------------------------
            class_source = get_most_recent_file(rec_dir, 'class_*.txt');
            if ~isempty(class_source)
                class = load_classification(class_source);
                fprintf('  %s: Loaded classification from %s\n', datestr(now), class_source);
                assert(length(class)==obj.num_cells,...
                       sprintf('Number of labels in %s is not consistent with %s!',...
                               class_source, data_source));
            else % No classification file
                fprintf('  %s: No classification file in %s!\n', datestr(now), rec_dir);
                class = cell(obj.num_cells,1); % Empty
            end

            images = squeeze(num2cell(data.filters, [1 2])); % images{k} is the 2D image of cell k
            [height, width] = size(images{1});
            boundaries = cell(size(images));
            masks = cell(size(images));
            for k = 1:obj.num_cells
                boundary = compute_ic_boundary(images{k}, 0.3);
                boundaries{k} = boundary{1}; % Keeps only the longest-boundary!
                masks{k} = poly2mask(boundaries{k}(:,1), boundaries{k}(:,2),...
                                     height, width);
            end
            
            obj.cells = struct(...
                'im', images,...
                'boundary', boundaries,...
                'mask', masks,...
                'label', class);
        end
        
        % Helper functions
        %------------------------------------------------------------
        function turn = compute_turn(obj, start, final)
            path = {start, final};
            if (all(strcmp(path, {'east', 'south'})) || ...
                all(strcmp(path, {'south', 'west'})) || ...
                all(strcmp(path, {'west', 'north'})) || ...
                all(strcmp(path, {'north', 'east'})))
                turn = 'left';
            else
                turn = 'right';
            end
        end
        
        % Accessors
        %------------------------------------------------------------
        function [trace, frame_indices] = get_trace(obj, cell_idx, selected_trials)
            % When 'selected_trials' is omitted, then return all trials
            if ~exist('selected_trials', 'var')
                selected_trials = 1:obj.num_trials;
            end
            
            trace = [];
            frame_indices = [];
            for k = selected_trials
                trace = [trace obj.trials(k).traces(cell_idx,:)]; %#ok<*AGROW>
                frame_indices = [frame_indices obj.trial_indices(k,1):obj.trial_indices(k,end)];
            end
        end
        
        function mask = get_mask(obj, cell_indices)
            % When 'cell_indices' is omitted, then return the masks of all
            % classified cells
            if ~exist('cell_indices', 'var')
                cell_indices = find(obj.is_cell);
            end
            
            [height, width] = size(obj.cells(1).im);
            mask = zeros(height, width);
            for cell_idx = cell_indices
                mask = mask | obj.cells(cell_idx).mask;
            end
        end
        
        function is_cell = is_cell(obj, cell_indices)
            % When 'cell_indices' is omitted, then return the label of all
            % cells
            if ~exist('cell_indices', 'var')
                cell_indices = 1:obj.num_cells;
            end
            
            is_cell = zeros(size(cell_indices));
            for k = 1:length(cell_indices)
                cell_idx = cell_indices(k);
                is_cell(k) = any(strcmp(obj.cells(cell_idx).label,...
                    {'phase-sensitive cell', 'cell'}));
            end
        end
        
        function is_correct = get_trial_correctness(obj)
            is_correct = cellfun(@strcmp, {obj.trials.goal}, {obj.trials.end});
        end
        
        % Debug functions
        %------------------------------------------------------------
        function set_all_labels(obj, label)
            % Overwrites the label of all cells in DaySummary
            for k = 1:obj.num_cells
                obj.cells(k).label = label;
            end
        end
        
        % Built-in visualization functions
        % Note: Do NOT make use of subplots in the built-in plot methods
        %------------------------------------------------------------
        function plot_cell_map(obj, color_grouping)
            % Optional argument allows for specification of color used for
            % the cell in the cell map. The color specification is defined
            % as follows:
            %   color_grouping = {[1, 2, 3, 4], [5, 6], [10]}
            % means that cells [1, 2, 3, 4] will be displayed in one color,
            % cells [5, 6] in another color, and [10] in another.
            
            % By default, color the cells based on classification
            if ~exist('color_grouping', 'var')
                cell_colors = arrayfun(@num2color, obj.is_cell());
            else
                cell_colors = repmat('w', 1, obj.num_cells); % Ungrouped cells are white
                % Unpack the colors
                for k = 1:length(color_grouping)
                    for cell_idx = color_grouping{k}
                        cell_colors(cell_idx) = num2color(k);
                    end
                end
            end
            
            % Background image to display
            [height, width] = size(obj.cells(1).im);
            ref_image = zeros(height, width);
            for k = 1:obj.num_cells
                ref_image = ref_image + obj.cells(k).im;
            end
            imagesc(ref_image);
            colormap gray;
            axis equal;
            xlim([1 width]);
            ylim([1 height]);
            
            hold on;
            for k = 1:obj.num_cells
                color = cell_colors(k);
                boundary = obj.cells(k).boundary;

                plot(boundary(:,1), boundary(:,2), 'Color', color);
%                 text(max(boundary(:,1)), min(boundary(:,2)),...
%                      sprintf('%d', k), 'Color', color);
            end
            hold off;
            
            function color = num2color(num)
                switch num
                    case 0
                        color = 'r';
                    case 1
                        color = 'g';
                    case 2
                        color = 'm';
                end
            end
        end
        
        function plot_trace(obj, cell_idx)
            % Plot the trace of a single cell, color-coded by trial
            trace_min = Inf;
            trace_max = -Inf;
            
            colors = 'kbr';
            for k = 1:obj.num_trials
                [trace, inds] = obj.get_trace(cell_idx, k);
                
                trial_trace_min = min(trace);
                trial_trace_max = max(trace);
                if (trial_trace_min < trace_min)
                    trace_min = trial_trace_min;
                end
                if (trial_trace_max > trace_max)
                    trace_max = trial_trace_max;
                end
                
                plot(inds, trace,...
                     colors(mod(k,length(colors))+1));
                hold on;
            end
            hold off;
            xlim([1 inds(end)]); % Using the last 'inds' from loop!
            trace_delta = trace_max - trace_min;
            ylim([trace_min trace_max] + 0.1*trace_delta*[-1 1]);
            xlabel('Frame index');
            ylabel('Signal [a.u.]');
        end
        
        function filtered_trials = filter_trials(obj, varargin)
            filtered_trials = ones(1, obj.num_trials);
            for k = 1:length(varargin)
                vararg = varargin{k};
                if ischar(vararg)
                    switch lower(vararg)
                        case 'incorrect'
                            filtered_trials = filtered_trials &...
                                (~strcmp({obj.trials.goal}, {obj.trials.end}));
                        case 'correct'
                            filtered_trials = filtered_trials &...
                                strcmp({obj.trials.goal}, {obj.trials.end});
                        case 'start'
                            filtered_trials = filtered_trials &...
                                strcmp({obj.trials.start}, varargin{k+1});
                        case 'end'
                            filtered_trials = filtered_trials &...
                                strcmp({obj.trials.end}, varargin{k+1});
                        case 'turn'
                            filtered_trials = filtered_trials &...
                                strcmp({obj.trials.turn}, varargin{k+1});
                    end
                end
            end
            filtered_trials = filtered_trials';
        end
        
        function plot_superposed_trials(obj, cell_idx, varargin)
            % Optional arguments allow for filtering of trials, e.g.
            %   "plot_superposed_trials(cell_idx, 'start', 'east')"
            display_trial = ones(obj.num_trials, 1);
            if ~isempty(varargin)
                display_trial = obj.filter_trials(varargin{:});
            end
            
            trace_min = Inf;
            trace_max = -Inf;
            
            colors = 'kbr';
            for k = 1:obj.num_trials
                if display_trial(k)
                    trial_trace = obj.get_trace(cell_idx, k);

                    trial_trace_min = min(trial_trace);
                    trial_trace_max = max(trial_trace);
                    if (trial_trace_min < trace_min)
                        trace_min = trial_trace_min;
                    end
                    if (trial_trace_max > trace_max)
                        trace_max = trial_trace_max;
                    end

                    plot(linspace(0, 1, length(trial_trace)),...
                         trial_trace,...
                         colors(mod(k,length(colors))+1));
                    hold on;
                end
            end
            hold off;
            grid on;
            xlim([0 1]);
            trace_delta = trace_max - trace_min;
            ylim([trace_min trace_max] + 0.1*trace_delta*[-1 1]);
            xlabel('Trial phase [a.u.]');
            ylabel('Trace [a.u.]');
        end
        
        function raster = plot_cell_raster(obj, cell_idx, varargin)
            % Optional argument 'draw_correct' will place a box at the end
            % of each trial indicating correct (green) or incorrect (red)
            %
            % Additional optional arguments allow for filtering of trials,
            % e.g. "plot_cell_raster(cell_idx, 'start', 'east')"
            
            display_trial = ones(obj.num_trials, 1);
            draw_correct = 0;
            if ~isempty(varargin)
                for k = 1:length(varargin)
                    vararg = varargin{k};
                    if ischar(vararg)
                        switch lower(vararg)
                            case 'draw_correct'
                                draw_correct = 1;
                        end
                    end
                end
                % Trial filtering arguments
                display_trial = obj.filter_trials(varargin{:});
            end
                        
            resample_grid = linspace(0, 1, 1000);
            num_filtered_trials = sum(display_trial);
            raster = zeros(num_filtered_trials, length(resample_grid));
            correctness = zeros(num_filtered_trials);
            counter = 0;
            for k = 1:obj.num_trials
                if display_trial(k)
                    counter = counter+1;
                    line = obj.get_trace(cell_idx, k);
                    raster(counter,:) = interp1(linspace(0,1,length(line)),...
                                      line,...
                                      resample_grid,...
                                      'pchip');
                    correctness(counter) = obj.trials(k).correct;
                end
            end
            
            imagesc(resample_grid, 1:size(raster,1), raster);
            colormap jet;
            xlabel('Trial phase [a.u.]');
            xlim([0 1]);
            ylabel('Trial index');
            
            if (draw_correct)
                corr_width = 0.025;
                for k = 1:num_filtered_trials
                    if correctness(k)
                        corr_color = 'g';
                    else
                        corr_color = 'r';
                    end
                    rectangle('Position', [1 k-0.5 corr_width 1],...
                              'FaceColor', corr_color);
                end
                xlim([0 1+corr_width]);
            end
        end
    end
end

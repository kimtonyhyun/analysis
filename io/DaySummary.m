% Summary of PlusMaze data for a single day.
%
% Inputs:
%   plusmaze_txt: Name of the PlusMaze output text file. Frame indices of
%       the text file needs to be consistent with the ICA traces!
%   ica_dir: Directory containing 
%       - ICA results in a "ica_*.mat" file
%       - Classification results in a "class_*.txt" file
%
% Output:
%   DaySummary object
%
% Example usage:
%   ds = DaySummary('c9m7d06_ti2.txt', 'ica001');
%
classdef DaySummary
    properties
        cells
        trials
        
        num_cells
        num_trials
        
        trial_indices
    end
    
    methods
        function obj = DaySummary(plusmaze_txt, ica_dir, varargin)
            % Handle optional input
            exclude_probe_trials = 0;
            use_reconstruction = 0;
            for k = 1:length(varargin)
                if ischar(varargin{k})
                    switch lower(varargin{k})
                        case 'excludeprobe'
                            exclude_probe_trials = 1;
                        case 'reconst'
                            use_reconstruction = 1;
                    end
                end
            end
            
            % Load data
            %------------------------------------------------------------
            [trial_indices, loc_info, trial_durations] =...
                parse_plusmaze(plusmaze_txt);
            
            if use_reconstruction
                data_source = get_most_recent_file(ica_dir, 'rec_*.mat');
            else
                data_source = get_most_recent_file(ica_dir, 'ica_*.mat');
            end
            data = load(data_source);
            obj.num_cells = data.info.num_pairs;
            fprintf('%s: Loaded data from %s\n', datestr(now), data_source);
            
            % Parse trial data
            %   TODO: Bring in centroids corresponding to mouse position
            %------------------------------------------------------------
            if (exclude_probe_trials)
                is_probe = strcmp(loc_info(:,1), 'north') | ...
                           strcmp(loc_info(:,1), 'south');
                       
                trial_indices = trial_indices(~is_probe,:);
                loc_info = loc_info(~is_probe,:);
                trial_durations = trial_durations(~is_probe);
            end
            
            num_trials = size(trial_indices, 1);
            traces = cell(num_trials, 1);
            for k = 1:num_trials
                trial_frames = trial_indices(k,1):...
                               trial_indices(k,end);
                traces{k} = data.traces(trial_frames, :)';
            end
            
            obj.num_trials = num_trials;
            obj.trial_indices = trial_indices;
            obj.trials = struct(...
                'start', loc_info(:,1),...
                'goal',  loc_info(:,2),...
                'end',   loc_info(:,3),...
                'time',  num2cell(trial_durations),...
                'traces', traces);
            
            % Parse cell data
            %   TODO: Consider classification as an optional input
            %------------------------------------------------------------
            class_source = get_most_recent_file(ica_dir, 'class_*.txt');
            class = load_classification(class_source);
            assert(length(class)==obj.num_cells,...
                   sprintf('Number of labels in %s is not consistent with %s!',...
                           class_source, data_source));

            images = squeeze(num2cell(data.filters, [1 2])); % images{k} is the 2D image of cell k
            boundaries = cell(size(images));
            for k = 1:obj.num_cells
                boundary = compute_ic_boundary(images{k}, 0.3);
                boundaries{k} = boundary{1};
            end
            
            obj.cells = struct(...
                'im', images,...
                'boundary', boundaries,...
                'label', class);
        end
        
        function [trace, frame_indices] = get_trace(obj, cell_idx, varargin)
            % Optional varargin specifies subset of trials. If omitted,
            % then pull the trace from all trials
            if isempty(varargin)
                selected_trials = 1:obj.num_trials;
            else
                selected_trials = varargin{1};
            end
            
            trace = [];
            frame_indices = [];
            for k = selected_trials
                trace = [trace obj.trials(k).traces(cell_idx,:)]; %#ok<*AGROW>
                frame_indices = [frame_indices obj.trial_indices(k,1):obj.trial_indices(k,end)];
            end
        end
        
        function is_cell = is_cell(obj, varargin)
            % Optional argument specifies subset of cell indices to check
            cell_indices = 1:obj.num_cells;
            if ~isempty(varargin)
                cell_indices = varargin{1};
            end
            is_cell = ~strcmp({obj.cells(cell_indices).label}, 'not a cell');
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
                if obj.is_cell(k)
                    color = 'g';
                else
                    color = 'r';
                end
                plot(boundary(:,1), boundary(:,2), 'Color', color);
                text(max(boundary(:,1)), min(boundary(:,2)),...
                     sprintf('%d', k), 'Color', color);
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
            % Optional arguments allow for filtering of trials, e.g.
            %   "plot_cell_raster(cell_idx, 'start', 'east')"
            display_trial = ones(obj.num_trials, 1);
            if ~isempty(varargin)
                display_trial = obj.filter_trials(varargin{:});
            end
                        
            resample_grid = linspace(0, 1, 1000);
            raster = zeros(sum(display_trial), length(resample_grid));            
            counter = 0;
            for k = 1:obj.num_trials
                if display_trial(k)
                    counter = counter+1;
                    line = obj.get_trace(cell_idx, k);
                    raster(counter,:) = interp1(linspace(0,1,length(line)),...
                                      line,...
                                      resample_grid,...
                                      'pchip');
                end
            end
            
            imagesc(resample_grid, 1:size(raster,1), raster);
            xlabel('Trial phase [a.u.]');
            ylabel('Trial index');
        end
    end
end

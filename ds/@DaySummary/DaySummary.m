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
    
    properties (SetAccess = private, Hidden)
        cell_map_ref_img
    end
        
    methods
        function obj = DaySummary(plusmaze_txt, rec_dir, varargin)
            % Handle optional input
            exclude_probe_trials = 0;
            for k = 1:length(varargin)
                if ischar(varargin{k})
                    switch lower(varargin{k})
                        case {'excludeprobe', 'noprobe'}
                            exclude_probe_trials = 1;
                    end
                end
            end
            
            % Parse trial data
            %------------------------------------------------------------
            [trial_indices, loc_info, trial_durations] =...
                parse_plusmaze(plusmaze_txt); %#ok<*PROP>
            fprintf('  %s: Loaded trial metadata from %s\n', datestr(now), plusmaze_txt);
            
            % Load data
            %------------------------------------------------------------
            data_source = get_most_recent_file(rec_dir, 'rec_*.mat');
            data = load(data_source);
            obj.num_cells = data.info.num_pairs;
            fprintf('  %s: Loaded filters and traces from %s\n', datestr(now), data_source);
            
            % Check that the length of traces is consistent with the table
            % of trial indices.
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
            class = cell(obj.num_cells,1);
            
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
            
            % Load classification
            %------------------------------------------------------------
            class_source = get_most_recent_file(rec_dir, 'class_*.txt');
            if ~isempty(class_source)
                obj.load_class(class_source);
                fprintf('  %s: Loaded classification from %s\n', datestr(now), class_source);
            end
            
            % Precompute cell map image, to avoid doing it each time
            %------------------------------------------------------------
            [height, width] = size(obj.cells(1).im);
            ref_image = zeros(height, width);
            for k = 1:obj.num_cells
                ref_image = ref_image + obj.cells(k).im;
            end
            obj.cell_map_ref_img = ref_image;
        end
        
        % Helper functions
        %------------------------------------------------------------
        function turn = compute_turn(~, start, final)
            % TODO: Turn into Static
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
        
        % Classification
        %------------------------------------------------------------
        function class = get_class(obj)
            class = {obj.cells.label}'; % Column cell
        end
        
        function set_all_labels(obj, label)
            % Overwrites the label of all cells in DaySummary
            for k = 1:obj.num_cells
                obj.cells(k).label = label;
            end
        end
    end
end

% Summary of PlusMaze data for a single day.
%
% Inputs:
%   ds_source: Three possibilities:
%       1) Plus maze text file
%       2) Struct 's' containing the following fields:
%           - s.maze: Path to plus maze text file (required)
%           - s.behavior: Path to behavioral video (optional; e.g. mp4)
%           - s.tracking: Path to tracking text file (optional; *.xy)
%       3) Empty string. Fake trial metadata will be generated that matches
%           the number of frames in the provided rec file. A hack to allow
%           use of the 'analysis' codebase with non-PlusMaze datasets.
%
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
        switchdata
        
        num_cells
        num_trials
        
        trial_indices
        full_num_frames
    end
    
    properties (SetAccess = private, Hidden=true)
        trace_corrs
        cell_distances
        cell_map_ref_img
        behavior_vid
        behavior_ref_img
        is_tracking_loaded
        is_eventdata_loaded
        orig_trial_indices
    end
        
    methods
        function obj = DaySummary(ds_source, rec_dir, varargin)
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
            
            % Load extraction data (i.e. filters & traces)
            %------------------------------------------------------------
            data_source = get_most_recent_file(rec_dir, 'rec_*.mat');
            data = load(data_source);
            obj.num_cells = data.info.num_pairs;
            fprintf('%s: Loaded %d filters and traces (%s) from %s\n',...
                    datestr(now), obj.num_cells, data.info.type, data_source);
                
            trace_num_frames = size(data.traces, 1);
            
            % Check if DaySummary session data is provided as a struct
            %------------------------------------------------------------
            if isstruct(ds_source)
                plusmaze_txt = ds_source.maze;
            else
                plusmaze_txt = ds_source;
            end
                       
            if ~isempty(plusmaze_txt)
                % PlusMaze metadata provided. Read trial metadata
                [trial_indices, loc_info, trial_durations] = parse_plusmaze(plusmaze_txt); %#ok<*PROP>
                fprintf('%s: Loaded trial metadata from %s\n', datestr(now), plusmaze_txt);
            else
                % PlusMaze metadata NOT provided. Make up fake info so that
                % DaySummary object can be instantiated anyway.
                trial_indices = [1 2 3 trace_num_frames];
                loc_info = {'east', 'north', 'north'};
                trial_durations = 10.0;
                fprintf('%s: Generated fake trial metadata to match number of frames (%d) in rec file (%s)!\n',...
                    datestr(now), trace_num_frames, data_source);
            end
            
            % Check that the length of traces is consistent with the table
            % of trial indices.
            obj.full_num_frames = trial_indices(end,end);
            assert(trace_num_frames == obj.full_num_frames,...
                'Error: Length of traces does not match trial index table!');
            
            % Optional exclusion of probe trials. Effectively, we are
            % "deleting" the lines of the plus maze text file that
            % correspond to probe trials.
            if (exclude_probe_trials)
                is_probe = strcmp(loc_info(:,1), 'north') | ...
                           strcmp(loc_info(:,1), 'south');
                       
                trial_indices = trial_indices(~is_probe,:);
                loc_info = loc_info(~is_probe,:);
                trial_durations = trial_durations(~is_probe);
            end
            
            % Parse by TRIAL
            %------------------------------------------------------------
            num_trials = size(trial_indices, 1);
            turns = cell(num_trials, 1);
            traces = cell(num_trials, 1);
            centroids = cell(num_trials, 1);

            for k = 1:num_trials
                trial_frames = trial_indices(k,1):...
                               trial_indices(k,end);

                num_frames_in_trial = length(trial_frames);
                traces{k} = data.traces(trial_frames, :)';
                turns{k} = obj.compute_turn(loc_info{k,1}, loc_info{k,3});
                centroids{k} = zeros(num_frames_in_trial, 2);
            end
            
            obj.num_trials = num_trials;
            
            % NOTE: We have to be quite careful about how we handle trial
            % indices if we apply probe trial elimination.
            %
            % Basically, when we need to refer to the raw data, we need the
            % _original_ frame indices, since the raw data sources do not
            % omit probe trials.
            %
            % On the other hand, when we later access content from the
            % probe-trial removed DaySummary (e.g. get_trace), the
            % resulting content does not know about the omitted trials.
            % Hence, for self-consistency of the DaySummary instance, we
            % need to "compress" the trial indices to remove gaps due to 
            % the omitted probe trials.
            %
            obj.trial_indices = compress_frame_indices(trial_indices, [0 0]);
            obj.orig_trial_indices = double(trial_indices);
            obj.trials = struct(...
                'start', loc_info(:,1),...
                'goal',  loc_info(:,2),...
                'end',   loc_info(:,3),...
                'correct', cellfun(@strcmp, loc_info(:,2), loc_info(:,3), 'UniformOutput', false),...
                'turn',  turns,...
                'time',  num2cell(trial_durations),...
                'traces', traces,...
                'events', [],...
                'centroids', centroids);
            
            % Compute trace correlation among all sources
            full_traces = cell2mat(traces'); % [num_cells x num_frames]
            fprintf('  Computing trace correlations between all sources...');
            tic;
            obj.trace_corrs = corr(full_traces');
            t = toc;
            fprintf(' Done (%.1f sec)\n', t);
            
            % Parse by CELL
            %------------------------------------------------------------
            class = cell(obj.num_cells,1);
            
            images = squeeze(num2cell(data.filters, [1 2])); % images{k} is the 2D image of cell k
            boundaries = cell(size(images));
            masks = cell(size(images));
            coms = cell(size(images)); % Center of mass

            fprintf('  Computing auxiliary spatial parameters...');
            tic;
            [height, width] = size(images{1});
            for k = 1:obj.num_cells
                boundary = compute_ic_boundary(images{k}, 0.3);
                boundaries{k} = boundary{1}; % Keep only the longest boundary!
                masks{k} = poly2mask(boundaries{k}(:,1), boundaries{k}(:,2), height, width);
                
                % Compute the center of mass
%                 masked_filter = masks{k}.*images{k};
                masked_filter = images{k}; % Masking sometimes yields numerical instabilities
                com = [(1:width)*sum(masked_filter,1)';
                       (1:height)*sum(masked_filter,2)];
                coms{k} = com / sum(masked_filter(:));
            end
            t = toc;
            fprintf(' Done (%.1f sec)\n', t);
            
            obj.cells = struct(...
                'im', images,...
                'boundary', boundaries,...
                'mask', masks,...
                'com', coms,...
                'label', class);
            
            % Compute distances among all sources
            fprintf('  Computing distances between all sources...');
            tic;
            D = Inf*ones(obj.num_cells);
            for i = 1:(obj.num_cells-1)
                for j = (i+1):obj.num_cells
                    delta = obj.cells(i).com - obj.cells(j).com;
                    D(i,j) = norm(delta);
                end
            end
            obj.cell_distances = min(D, D'); % Make symmetric. Note Infs.
            t = toc;
            fprintf(' Done (%.1f sec)\n', t);
            
            % Precompute cell map image, to avoid doing it each time
            [height, width] = size(obj.cells(1).im);
            ref_image = zeros(height, width);
            for k = 1:obj.num_cells
                ref_image = ref_image + obj.cells(k).im;
            end
            obj.cell_map_ref_img = ref_image;
            
            % Load classification, if available
            %------------------------------------------------------------
            class_source = get_most_recent_file(rec_dir, 'class_*.txt');
            if ~isempty(class_source)
                obj.load_class(class_source);
                fprintf('%s: Loaded classification from %s\n', datestr(now), class_source);
            end
                       
            % Other initialization
            %------------------------------------------------------------
            obj.behavior_vid = [];
            obj.is_tracking_loaded = false;
            
            % Load behavior video and tracking data if available
            if isstruct(ds_source)
                if isfield(ds_source, 'behavior')
                    obj.load_behavior_movie(ds_source.behavior);
                end
                if isfield(ds_source, 'tracking')
                    obj.load_tracking(ds_source.tracking);
                end
            end
            
            % Event data
            obj.is_eventdata_loaded = false;
            event_source = get_most_recent_file(rec_dir, 'events_*.mat');
            if ~isempty(event_source)
                obj.load_events(event_source);
            end
            
            % Fill in switch data
            obj.switchdata = struct(...
                'pre_switch_trials', [],...
                'post_switch_trials', [],...
                'changing_path_start', '',...
                'constant_path_start', '',...
                'changing_path_cutoff', []);
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
        end % compute_turn

        function filtered_trials = filter_trials(obj, varargin)
            filtered_trials = ones(1, obj.num_trials);
            for k = 1:length(varargin)
                vararg = varargin{k};
                if ischar(vararg)
                    switch lower(vararg)
                        case {'range', 'inds'}
                            % Convert list of trial indices to logical vector
                            lv = ismember(1:obj.num_trials, varargin{k+1});
                            filtered_trials = filtered_trials & lv;
                        case 'incorrect'
                            filtered_trials = filtered_trials &...
                                (~strcmp({obj.trials.goal}, {obj.trials.end}));
                        case 'correct'
                            filtered_trials = filtered_trials &...
                                strcmp({obj.trials.goal}, {obj.trials.end});
                        case 'start'
                            filtered_trials = filtered_trials &...
                                trial_filter(obj, {obj.trials.start}, varargin{k+1});
                        case 'end'
                            filtered_trials = filtered_trials &...
                                trial_filter(obj, {obj.trials.end}, varargin{k+1});
                        case 'turn'
                            filtered_trials = filtered_trials &...
                                trial_filter(obj, {obj.trials.turn}, varargin{k+1});
                    end
                end
            end
            filtered_trials = filtered_trials';

            % helper function to filter cell arrays of strings, i.e. enable
            % selectors such as: filter_trials('start', {'east', 'west'})
            function mask = trial_filter(~, trial_data, selection)
                if ischar(selection)
                    mask = strcmp(trial_data,selection);
                elseif iscell(selection)
                    mask = false(size(trial_data));
                    for a = 1:length(selection)
                        mask = mask | strcmp(trial_data,selection{a});
                    end
                else
                    error('selection criterion must be a cell array or string')
                end
            end

        end   

        % Accessors
        %------------------------------------------------------------
        function count = num_classified_cells(obj)
            count = sum(obj.is_cell);
        end
        
        function [trace, frame_indices, selected_trials] = get_trace(obj, cell_idx, varargin)
            selected_trials = 1:obj.num_trials;
            for k = 1:length(varargin)
                vararg = varargin{k};
                if isnumeric(vararg)
                    selected_trials = vararg;
                end
            end
            filtered_trials = find(obj.filter_trials(varargin{:}));
            selected_trials = intersect(selected_trials, filtered_trials)';
            
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
              
        function is_correct = get_trial_correctness(obj)
            is_correct = [obj.trials.correct];
        end
        
        function selected_idx = get_cell_by_xy(obj, xy, varargin)
            % Returns the first cell index whose filter boundary encloses
            % the XY position provided in 'xy'
            
            % By default all cell candidates can be clicked
            classified_cells_only = 0;
            for i = 1:length(varargin)
                if ischar(varargin{i})
                    switch lower(varargin{i})
                        case {'cells', 'cellsonly'}
                            classified_cells_only = 1;
                    end
                end
            end
            
            selected_idx = []; % If no hit, then return empty
            for k = 1:obj.num_cells
                boundary = obj.cells(k).boundary;
                if (obj.is_cell(k) || ~classified_cells_only)
                    if inpolygon(xy(1), xy(2), boundary(:,1), boundary(:,2))
                        selected_idx = k;
                        break;
                    end
                end
            end
        end
        
        function neighbor_inds = get_nearest_sources(obj, cell_idx, num_neighbors, varargin)
            % Return the indices of 'num_neighbors' number of sources
            % closest to the source with 'cell_idx'.
            
            classified_cells_only = 0;
            for i = 1:length(varargin)
                if ischar(varargin{i})
                    switch lower(varargin{i})
                        case {'cell', 'cells', 'cellsonly'}
                            classified_cells_only = 1;
                    end
                end
            end

            d = obj.cell_distances(cell_idx,:); % Distance to all cells
            if classified_cells_only
                is_not_cell = ~obj.is_cell;
                d(is_not_cell) = Inf; % Set distances to non-cells as Inf
            end
            [~, neighbor_inds] = sort(d); % Ascending order
            neighbor_inds = neighbor_inds(1:num_neighbors);
        end
        
        % Classification
        %------------------------------------------------------------
        function class = get_class(obj)
            class = {obj.cells.label}'; % Column cell
        end
        
        function apply_labels_to(obj, label, cell_indices)
            for k = cell_indices
                obj.cells(k).label = label;
            end
        end
        
        function set_labels(obj, varargin)
            cell_indices = 1:obj.num_cells;
            if ~isempty(varargin)
                cell_indices = varargin{1};
            end
            obj.apply_labels_to('cell', cell_indices);
        end
        
        function reset_labels(obj, varargin)
            cell_indices = 1:obj.num_cells;
            if ~isempty(varargin)
                cell_indices = varargin{1};
            end
            obj.apply_labels_to([], cell_indices);
        end
        
        function set_unlabeled_cells(obj)
            unlabeled_cells = find(cellfun(@isempty, obj.get_class))';
            obj.apply_labels_to('not a cell', unlabeled_cells);
        end
        
        function invert_labels(obj)
            orig_cells = find(obj.is_cell);
            orig_not_cells = setdiff(1:obj.num_cells, orig_cells);
            obj.apply_labels_to('not a cell', orig_cells);
            obj.apply_labels_to('cell', orig_not_cells);
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
        
        % Load behavior movie
        %------------------------------------------------------------
        function loaded = is_behavior_loaded(obj)
            loaded = ~isempty(obj.behavior_vid);
        end
        
        function load_behavior_movie(obj, behavior_source)
            obj.behavior_vid = VideoReader(behavior_source);
            fprintf('%s: Loaded behavior video from "%s"\n',...
                datestr(now), behavior_source);
            if (obj.behavior_vid.NumberOfFrames ~= obj.full_num_frames)
                fprintf('  Warning! Number of frames in behavior video (%d) does not match the trial frame table (%d)!\n',...
                    obj.behavior_vid.NumberOfFrames, obj.trial_indices(end,end));
            end
            
            % Load reference image
            img = obj.behavior_vid.read(1);
            obj.behavior_ref_img = squeeze(img(:,:,1,:));
        end
        
        function Mb = get_behavior_trial(obj, trial_idx)
            % Must use original frame indices to access raw data
            trial_frames = obj.orig_trial_indices(trial_idx, [1 end]); % [Start end]
            Mb = obj.behavior_vid.read(trial_frames);
            Mb = squeeze(Mb(:,:,1,:)); % Movie is actually grayscale!
        end
        
        % Load tracking data
        %------------------------------------------------------------
        function load_tracking(obj, tracking_source)
            centroids = load(tracking_source);

            % Must use original frame indices to access raw data
            for k = 1:obj.num_trials
                inds = obj.orig_trial_indices(k,1):obj.orig_trial_indices(k,end);
                obj.trials(k).centroids = centroids(inds, :);
            end
            
            fprintf('%s: Loaded tracking data from "%s"\n',...
                datestr(now), tracking_source);
            obj.is_tracking_loaded = true;
        end
        
        % Load event data
        %------------------------------------------------------------
        function load_events(obj, event_source)
            data = load(event_source);
            assert(length(data.events) == obj.num_cells,...
                'Error: Number of cells in event file does not match that in DaySummary!');
            
            assert(data.events(1).info.num_frames == obj.full_num_frames,...
                'Error: Number of frames in event file does not match full number of frames in DaySummary!');
            
            events_per_trial = compute_events_per_trial({data.events.auto}, obj.orig_trial_indices);
            for k = 1:obj.num_trials
                obj.trials(k).events = events_per_trial{k};
            end
            
            fprintf('%s: Loaded events from "%s"\n',...
                datestr(now), event_source);
            obj.is_eventdata_loaded = true;
        end
        
        % Inspect switch data
        %------------------------------------------------------------
        function loaded = is_switchdata_loaded(obj)
            valid_starts = {'east', 'west'};
            loaded = ~isempty(obj.switchdata.pre_switch_trials) &&...
                     ~isempty(obj.switchdata.post_switch_trials) &&...
                     any(strcmp(obj.switchdata.changing_path_start, valid_starts)) &&...
                     any(strcmp(obj.switchdata.constant_path_start, valid_starts));
        end
        
        function st = get_switch_trials(obj)
            % Provide trial indices for switch analysis. Note that we only
            % examine correct trials.
            %
            % NOTE: Consider setting once when switchdata is loaded
            st = struct('constant_pre', [],...
                        'constant_post', [],...
                        'changing_pre', [],...
                        'changing_post', []);
                    
            if obj.is_switchdata_loaded
                sd = obj.switchdata;
                st.constant_pre = find(obj.filter_trials(...
                    'range', sd.pre_switch_trials,...
                    'start', sd.constant_path_start,...
                    'correct'));
                st.constant_post = find(obj.filter_trials(...
                    'range', sd.post_switch_trials,...
                    'start', sd.constant_path_start,...
                    'correct'));
                st.changing_pre = find(obj.filter_trials(...
                    'range', sd.pre_switch_trials,...
                    'start', sd.changing_path_start,...
                    'correct'));
                st.changing_post = find(obj.filter_trials(...
                    'range', sd.post_switch_trials,...
                    'start', sd.changing_path_start,...
                    'correct'));
            end
        end
        
        function reset_switchdata(obj)
            obj.switchdata.pre_switch_trials = [];
            obj.switchdata.post_switch_trials = [];
            obj.switchdata.changing_path_start = '';
            obj.switchdata.constant_path_start = '';
            obj.switchdata.changing_path_cutoff = [];
        end
            
    end % public methods

end

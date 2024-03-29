classdef MultiDay < handle
    properties (SetAccess = private)
        valid_days % List of days in the MultiDay object
        valid_day_names
        num_days
        num_cells % Number of aligned cells through all days
        
        matched_indices
        sort_day % 'matched_indices' is sorted in ascending order according
                 % to cell IDs on 'sort_day'
    end
    
    properties (Access = private)
        ds
        match

        full_to_sparse
        source_id_offsets
        num_all_sources
    end
    
    methods
        function obj = MultiDay(ds_list, match_list, varargin)
            % By default, MD keeps only cells that show up on ALL days
            keepall = false;
            for k = 1:length(varargin)
                if ischar(varargin{k})
                    switch varargin{k}
                        % "Keep-all" option:
                        % Keeps cells if it shows up on at least one
                        % DaySummary. USE WITH CAUTION!
                        case 'keepall'
                            keepall = true;
                            
                    end
                end
            end
            
            % Unpack the provided list of DaySummary's into a full cell
            % TODO: Consider using sparse internal storage.
            %------------------------------------------------------------
            obj.valid_days = cell2mat(ds_list(:,1))';
            obj.num_days = length(obj.valid_days);
            
            % Manually provided description of each "day"
            obj.valid_day_names = cell(1, obj.num_days);
            if size(ds_list, 2) >= 3 % ds_list has at least 3 columns
                for k = 1:obj.num_days
                    obj.valid_day_names{k} = ds_list{k,3};
                end
            else
                for k = 1:obj.num_days
                    obj.valid_day_names{k} = sprintf('Day %d', obj.valid_days(k));
                end
            end
            
            max_day = max(obj.valid_days);
            
            obj.full_to_sparse = zeros(1, max_day);
            obj.ds = cell(max_day, 1);
            for k = 1:size(ds_list,1)
                day  = ds_list{k,1};
                ds_k = ds_list{k,2};
                fprintf('%s: Day %d has %d classified cells (out of %d)\n',...
                    datestr(now), day, ds_k.num_classified_cells, ds_k.num_cells);
                
                obj.ds{day} = ds_k;
                obj.full_to_sparse(day) = k;
            end

            % Unpack the provided list of matches into a full cell matrix
            % Convention:
            %   match{i,j} = match from Day i to Day j
            %------------------------------------------------------------
            obj.match = cell(max_day, max_day);
            for k = 1:size(match_list,1)
                i = match_list{k,1};
                j = match_list{k,2};
                m_itoj = match_list{k,3};
                m_jtoi = match_list{k,4};
                
                if all(ismember([i, j], obj.valid_days))
                    obj.match{i,j} = m_itoj;
                    obj.match{j,i} = m_jtoi;
                end
            end
            
            % Check if match matrices have been provided
            for i = obj.valid_days
                for j = setdiff(obj.valid_days, i)
                    if isempty(obj.match{i,j})
                        fprintf('  Warning! No match provided from Day %d to Day %d!\n', i, j);
                    end
                end
            end
            
            % Compute matches. Basic idea is that each cell from each day
            % will be given a unique identifier (see 'day_to_linear_id'), 
            % and all pairwise match data will be incorporated into a
            % graph. We then look for connected components of the graph.
            %------------------------------------------------------------
            num_sources_per_day = zeros(1, obj.num_days);
            for k = 1:obj.num_days
                day = obj.valid_days(k);
                num_sources_per_day(k) = obj.ds{day}.num_cells; % Includes all sources, including non-cells
            end
            obj.num_all_sources = sum(num_sources_per_day);
            obj.source_id_offsets = cumsum([0 num_sources_per_day(1:end-1)]);

            A = obj.compute_adjacency_matrix;

            M = obj.compute_all_matches(A);
            if ~keepall
                M = obj.filter_matches(M); % Cell must be matched on every day
                M = obj.verify_match_list_consistency(M);
            end
            
            obj.matched_indices = M;
            obj.sort_matches_by_day(obj.valid_days(1));
            obj.num_cells = size(obj.matched_indices, 1);
            fprintf('%s: Found %d matching classified cells across all days\n',...
                datestr(now), obj.num_cells);
            
        end % MultiDay
        
        % Generate a new MultiDay object from a subset of days
        function nobj = get_sub_multiday(obj, sub_days)
            % Verify that all days are contained in the parent md
            if ~all(ismember(sub_days, obj.valid_days))
                error('Error! Parent MultiDay does not have the requested subset of days');
            end
            
            num_sub_days = length(sub_days);
            ds_list = cell(num_sub_days, 2);
            for i = 1:num_sub_days
                day = sub_days(i);
                ds_list(i,:) = {day, obj.ds{day}};
            end
            
            num_matches = nchoosek(num_sub_days,2);
            match_list = cell(num_matches, 4);
            idx = 1;
            for i = 1:(num_sub_days-1)
                day_i = sub_days(i);
                for j = (i+1):num_sub_days
                    day_j = sub_days(j);
                    match_list(idx,:) = {day_i, day_j, obj.match{day_i, day_j}, obj.match{day_j, day_i}};
                    idx = idx + 1;
                end
            end
            
            nobj = MultiDay(ds_list, match_list);
        end
        
        % General accessors
        %------------------------------------------------------------
        function ds = day(obj, day_idx)
            if ~ismember(day_idx, obj.valid_days)
                error('Error! Day %d is not valid for this MultiDay', day_idx);
            end
            ds = obj.ds{day_idx};
        end
        
        function indices = get_indices(obj, selected_days)
            % Returned indices are day-specific
            if ~exist('selected_days', 'var')
                selected_days = obj.valid_days;
            end
            selected_days = obj.full_to_sparse(selected_days);
            indices = obj.matched_indices(:, selected_days);
        end
        
        function trials = get_trials(obj, day_idx, varargin) 
            
            % Defaults
            ds = obj.day(day_idx); %#ok<*PROP>
            num_trials = ds.num_trials;
            trial_indices = (1:num_trials)';% get all trials
            filtered_trials = true(num_trials,1);
            split_trial_phases = 0;
            
            if ~isempty(varargin)
                for k = 1:length(varargin)
                    if ischar(varargin{k})
                        switch varargin{k}
                            case 'trials'
                                if isnumeric(varargin{k+1})
                                    trial_indices = varargin{k+1};
                                    if size(trial_indices,1)==1 % row vector
                                        trial_indices = trial_indices';
                                    end
                                end
                            case {'start', 'end', 'turn'}
                                var1 = varargin{k};
                                var2 = varargin{k+1};
                                if ~isempty(var2)
                                    filtered_trials = filtered_trials & ...
                                        ds.filter_trials(var1,var2);
                                end
                            case 'correct'
                                filtered_trials = filtered_trials & ...
                                    ds.filter_trials('correct');
                            case 'incorrect'
                                filtered_trials = filtered_trials & ...
                                    ds.filter_trials('incorrect');
                            case {'split_trial_phases'}
                                split_trial_phases = 1;
                        end
                    end
                end
            end
            
            trial_indices = intersect(trial_indices,find(filtered_trials));
            trials = ds.trials(trial_indices);  
            if split_trial_phases
                frame_indices = ds.trial_indices(trial_indices,:);
                % Remove offsets
                frame_indices = bsxfun(@minus,frame_indices,...
                    frame_indices(:,1)-1);
                frame_indices = frame_indices(trial_indices,:);
            end
            
            % Reorder the traces to match the common (matched) index
            day_cell_indices = obj.get_indices(day_idx);
            for k = 1:length(trials)
                trials(k).traces = trials(k).traces(day_cell_indices, :);
                if split_trial_phases
                    idx_pre = frame_indices(k,1):frame_indices(k,2);
                    idx_run = (frame_indices(k,2)+1):frame_indices(k,3);
                    idx_post = (frame_indices(k,3)+1):frame_indices(k,4);
                    x = trials(k).traces;
                    x = {x(:,idx_pre), x(:,idx_run), x(:,idx_post)};
                    trials(k).traces = x;
                end
            end
        end
        
        % Accessors using "common index" (i.e. indexing variable over the
        % matched cells -- not specific to day)
        %------------------------------------------------------------
        function cell_idx = get_cell_idx(obj, common_cell_idx, day_idx)
            % Convert the common_cell_idx into the day-specific index
            cell_idx = obj.matched_indices(common_cell_idx,...
                            obj.full_to_sparse(day_idx));
        end
        
        function cell = get_cell(obj, common_cell_idx, day_idx)
            cell_idx = obj.get_cell_idx(common_cell_idx, day_idx);
            cell = obj.ds{day_idx}.cells(cell_idx);
        end
        
        % Auxiliary methods
        %------------------------------------------------------------
        function sort_inds = sort_matches_by_day(obj, day_idx)
            sort_col = obj.full_to_sparse(day_idx);
            [mi, sort_inds] = sortrows(obj.matched_indices, sort_col);
            
            % Special handling of "0", which is a flag for a lack of match.
            % For certain applications, it is useful to stack these at the
            % end of the list.
            unmatched_rows = (mi(:,sort_col) == 0);
            matched_rows = ~unmatched_rows;
            
            obj.matched_indices = [mi(matched_rows,:);
                                   mi(unmatched_rows,:)];
            sort_inds = [sort_inds(matched_rows);
                         sort_inds(unmatched_rows)];
            obj.sort_day = day_idx;
        end
        
        function unmatched_ids = get_unmatched_cells(obj, day_idx)
            % Retrieve the IDs of cells on 'day_idx' that are not matched
            % across all days of MultiDay. NOTE: The definition of
            % "unmatched" computed here considers cells that match to 
            % non-cells as unmatched!
            all_cell_ids = find(obj.day(day_idx).is_cell);
            matched_ids = obj.matched_indices(:, obj.full_to_sparse(day_idx));
            unmatched_ids = setdiff(all_cell_ids, matched_ids);
        end

        function trial_map = filter_trials(md,varargin)
            % TRIAL_MAP = FILTER_TRIALS(MD, OPTS...)
            %
            % Returns a num_trials x 2 matrix of integers specifying day number and
            % trial number for trials that meet specified criteria. For example, the
            % following returns trials starting in east arm, ending in north arm:
            %
            % trial_map = filter_trials(md, 'start', 'east', 'end', 'north')
            
            % preallocate
            max_trials = 0;
            for d = md.valid_days
                max_trials = max_trials + length(md.day(d).trials);
            end
            trial_map = zeros(max_trials,2);
            
            % filter each day
            a = 1;
            for d = md.valid_days
                trial_idx = find(md.day(d).filter_trials(varargin{:}));
                b = a-1 + length(trial_idx);
                trial_map(a:b,1) = d;
                trial_map(a:b,2) = trial_idx;
                a = b+1;
            end
            
            % truncate unused storage
            trial_map = trial_map(1:b,:);
        end

        % Functions for computing matches
        %------------------------------------------------------------
        function A = compute_adjacency_matrix(obj)
            A = false(obj.num_all_sources, obj.num_all_sources);
            for day_i = obj.valid_days
                for day_j = setdiff(obj.valid_days, day_i)
                    m_itoj = obj.match{day_i, day_j};
                    if ~isempty(m_itoj)
                        for cell_i = 1:length(m_itoj)
                            m = m_itoj{cell_i};
                            if ~isempty(m)
                                cell_j = m(1);
                                I = obj.day_to_linear_id(day_i, cell_i);
                                J = obj.day_to_linear_id(day_j, cell_j);
                                A(I,J) = true;
                            end
                        end
                    end
                end % day_j
            end % day_i

            if ~issymmetric(A)
                cprintf('Warning: Adjacency matrix is not symmetric!\n');
            end
        end

        function M = compute_all_matches(obj, A)
            G = graph(A); % Generate graph object from adjacency matrix
            ccs = conncomp(G, 'OutputForm', 'cell');
            num_ccs = length(ccs);

            % Set up the matrix of matched indices. Note: zero indicates
            % unmatched (or not classified to be a cell)
            M = zeros(num_ccs, obj.num_days);
            for i = 1:num_ccs
                cc = ccs{i}; % i-th connected component
                for linear_id = cc % linear id within cc
                    [k, cell_idx] = obj.linear_to_day_id(linear_id);
                    day = obj.valid_days(k);
                    if obj.ds{day}.is_cell(cell_idx)
                        if (M(i,k) ~= 0) % There is an existing assignment!
                            fprintf('  Warning! Cells %d and %d from Day %d match to the same cross-day set!\n',...
                                M(i,k), cell_idx, day);
                        end
                        M(i,k) = cell_idx;
                    end
                end
            end
        end % compute_all_matches

        function Mf = filter_matches(~, M)
            % Filter out rows of the match matrix M that has zeros
            % (indicating non-matched cells)
            unmatched = any(M==0, 2);
            Mf = M(~unmatched,:);
        end

        function Mf = verify_match_list_consistency(obj, M)
            % For each row of the match matrix M, explicitly check that the
            % result is consistent with every pairwise match data that was
            % originally provided.
            num_matches = size(M,1);
            is_valid = true(num_matches, 1);
            
            day_inds = 1:size(M,2);
            for di = day_inds % Source day
                day_i = obj.valid_days(di);
                for dj = setdiff(day_inds, di) % Target day
                    day_j = obj.valid_days(dj);
                    m_itoj = obj.match{day_i, day_j};
                    if ~isempty(m_itoj) % Day i to Day j match was provided
                        for k = 1:num_matches
                            cell_i = M(k,di);
                            cell_j = M(k,dj);
                            if (cell_i ~= 0) && (cell_j ~= 0) % There is a pairing in M
                                m = m_itoj{cell_i};
                                if isempty(m) % No corresponding match in m_itoj
                                    is_valid(k) = false;
                                else
                                    if (m(1) ~= cell_j) % The match in m_itoj disagrees
                                        is_valid(k) = false;
                                    end
                                end
                            end
                        end % k
                    end
                end % dj
            end % di
            
            num_invalid_matches = sum(~is_valid);
            Mf = M(is_valid, :);
            if num_invalid_matches > 0
                cprintf('red', '  Removed %d matches inconsistent with provided match_list\n', num_invalid_matches);
            end
        end % verify_match_consistency

        function linear_id = day_to_linear_id(obj, day_idx, cell_idx)
            day_idx_sparse = obj.full_to_sparse(day_idx);
            linear_id = obj.source_id_offsets(day_idx_sparse) + cell_idx;
        end

        function [k, cell_idx] = linear_to_day_id(obj, linear_idx)
            k = find((linear_idx-obj.source_id_offsets)>0, 1, 'last');
            cell_idx = linear_idx - obj.source_id_offsets(k);
        end

    end % Public methods
end
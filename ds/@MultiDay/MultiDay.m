% Structure for accessing multiple aligned DaySummary's
%
classdef MultiDay < handle
    properties (SetAccess = private)
        valid_days % List of days in the MultiDay object
        num_days
        num_cells % Number of aligned cells through all days
        matched_indices
    end
    
    properties (Access = private)
        full_to_sparse
        ds
        match
    end
    
    methods
        function obj = MultiDay(ds_list, match_list)
            % Unpack the provided list of DaySummary's into a full cell
            % TODO: Consider using sparse internal storage.
            %------------------------------------------------------------
            obj.valid_days = cell2mat(ds_list(:,1))';
            obj.num_days = length(obj.valid_days);
            
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
                
                obj.match{i,j} = m_itoj;
                obj.match{j,i} = m_jtoi;
            end
            
            % Check if match matrices have been provided
            for i = obj.valid_days
                for j = setdiff(obj.valid_days, i)
                    if isempty(obj.match{i,j})
                        fprintf('  Warning! No match provided from Day %d to Day %d!\n', i, j);
                    end
                end
            end
            
            % Compute matches
            M = obj.compute_all_matches();
            
            % Filter out rows of M with unmatched indices (i.e. zeros) and
            % store result
            obj.matched_indices = obj.filter_matches(M);
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
        
        function trials = get_trials(obj, day_idx, trial_indices)
            if ~exist('trial_indices', 'var')
                trial_indices = 1:obj.day(day_idx).num_trials;
            end
            trials = obj.day(day_idx).trials(trial_indices);
            
            % Reorder the traces to match the common (matched) index
            day_cell_indices = obj.get_indices(day_idx);
            for k = 1:length(trials)
                trials(k).traces = trials(k).traces(day_cell_indices, :);
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
    end
    
    methods (Access=private)
        function [M, assignments] = compute_all_matches(obj)
            % Cells of all days will be arranged linearly for the graph
            % computation.
            num_cells_per_day = zeros(1, obj.num_days);
            for k = 1:obj.num_days
                day = obj.valid_days(k);
                num_cells_per_day(k) = obj.ds{day}.num_cells;
            end
            num_all_cells = sum(num_cells_per_day);
            offsets = cumsum([0 num_cells_per_day(1:(end-1))]);

            % Adjacency matrix with all cells arranged linearly
            A = zeros(num_all_cells, num_all_cells);
            for day_i = obj.valid_days
                for day_j = setdiff(obj.valid_days, day_i)
                    m_itoj = obj.match{day_i, day_j};
                    if ~isempty(m_itoj)
                        for cell_i = 1:length(m_itoj)
                            m = m_itoj{cell_i};
                            if ~isempty(m)
                                cell_j = m(1);
                                I = day_to_linear(day_i, cell_i);
                                J = day_to_linear(day_j, cell_j);
                                A(I,J) = 1;
                            end
                        end
                    end
                end % day_j
            end % day_i
            
            A = sparse(A); % Required by graphconncomp
            [num_components, assignments] = graphconncomp(A, 'Directed', 'false');
            
            % Set up the match matrix; zero indicates unmatched (or not
            % classified to be a cell)
            M = zeros(num_components, obj.num_days);
            for I = 1:length(assignments)
                assignment = assignments(I);
                [k, cell_idx] = linear_to_day(I);
                day = obj.valid_days(k);
                if obj.ds{day}.is_cell(cell_idx) % Keep only classified cells
                    if (M(assignment, k) ~= 0) % There is an existing assignment!
                        fprintf('  Warning! Cells %d and %d from Day %d match to the same cross-day set!\n',...
                            M(assignment, k), cell_idx, obj.valid_days(k));
                    end
                    M(assignment, k) = cell_idx;
                end
            end
            
            % Helper indexing functions
            function linear_idx = day_to_linear(day_idx, cell_idx)
                day_idx_sparse = obj.full_to_sparse(day_idx);
                linear_idx = offsets(day_idx_sparse) + cell_idx;
            end % day_to_linear
            
            function [k, cell_idx] = linear_to_day(linear_idx)
                k = find((linear_idx-offsets)>0, 1, 'last');
                cell_idx = linear_idx - offsets(k);
            end % linear_to_day

        end % compute_all_matches
        
        function Mf = filter_matches(~, M)
            % Filter out rows of the match matrix M that has zeros
            % (indicating non-matched cells)
            unmatched = any(M==0, 2);
            Mf = M(~unmatched,:);
        end % filter_matches
    end
end
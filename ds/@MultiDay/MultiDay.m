% Structure for accessing multiple aligned DaySummary's
%
classdef MultiDay < handle
    properties
        valid_days % List of days in the MultiDay object
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
            num_days = length(obj.valid_days);
            max_day = max(obj.valid_days);
            
            obj.ds = cell(max_day, 1);
            for k = 1:size(ds_list,1)
                day  = ds_list{k,1};
                ds_k = ds_list{k,2};
                fprintf('%s: Day %d has %d classified cells\n',...
                    datestr(now), day, sum(ds_k.is_cell));
                
                obj.ds{day} = ds_k;
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
            
            % Verify that there is is a match matrix for every pair of
            % provided DaySummarys
            for i = obj.valid_days
                for j = setdiff(obj.valid_days, i)
                    if isempty(obj.match{i,j})
                        error('No match between Day %d and Day %d provided!', i, j);
                    end
                end
            end
            
            % Compute matching across all provided days. We keep only the
            % _classified cells_ of each day
            %------------------------------------------------------------
            base_day = obj.valid_days(1);
            base_day_cells = find(obj.ds{base_day}.is_cell);
            base_day_num_cells = length(base_day_cells);
            
            % Scratch space (sparse) for matching indices
            obj.full_to_sparse = zeros(1, max_day);
            obj.full_to_sparse(base_day) = 1;
            
            M = zeros(base_day_num_cells, num_days);
            M(:,1) = base_day_cells;
            
            % Fill out M columnwise
            for k = 2:num_days
                prev_day = obj.valid_days(k-1);
                curr_day = obj.valid_days(k);
                obj.full_to_sparse(curr_day) = k;
                for x = 1:base_day_num_cells
                    % First, check that the row corresponds to a valid
                    % matched cell on the previous (k-1) day
                    prev_cell_idx = M(x,k-1);
                    if (prev_cell_idx ~= 0)
                        % Next, check if the cell from the previous (k-1)
                        % day matches to the current (k) day
                        m = obj.match{prev_day, curr_day}{prev_cell_idx};
                        if ~isempty(m)
                            curr_cell_idx = m(1);
                            % Finally, we save the matched result if it is
                            % a valid classified cell on the current day
                            if obj.ds{curr_day}.is_cell(curr_cell_idx)
                                M(x,k) = curr_cell_idx;
                            end
                        end
                    end
                end
            end
            
            % Filter out rows of M with unmatched indices (i.e. zeros) and
            % store result
            unmatched = any(M==0, 2);
            obj.matched_indices = M(~unmatched, :);
            obj.num_cells = size(obj.matched_indices, 1);
            fprintf('%s: Found %d matching classified cells across all days\n',...
                datestr(now), obj.num_cells);
            
        end % MultiDay
        
        % Accessors
        %------------------------------------------------------------
        function ds = day(obj, day_idx)
            if ~ismember(day_idx, obj.valid_days)
                error('Error! Day %d is not valid for this MultiDay', day_idx);
            end
            ds = obj.ds{day_idx};
        end
        
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
end
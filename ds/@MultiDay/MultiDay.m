% Structure for accessing multiple aligned DaySummary's
%
classdef MultiDay < handle
    properties
        days % List of days in the MultiDay object
        num_cells % Number of aligned cells through all days
        matched_indices
    end
    
    properties (Access = private)
        ds
        match
    end
    
    methods
        function obj = MultiDay(ds_list, match_list)
            % Unpack the provided list of DaySummarys into a full cell
            %------------------------------------------------------------
            obj.days = cell2mat(ds_list(:,1))';
            num_days = length(obj.days);
            max_day = max(obj.days);
            obj.ds = cell(max_day, 1);
            for k = 1:size(ds_list,1)
                day  = ds_list{k,1};
                ds_k = ds_list{k,2};
                fprintf('%s: Day %d ds has %d classified cells\n',...
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
            for i = obj.days
                for j = setdiff(obj.days, i)
                    if isempty(obj.match{i,j})
                        error('No match between Day %d and Day %d provided!', i, j);
                    end
                end
            end
            
            % Compute matching across all provided days. We keep only the
            % _classified cells_ of each day
            %------------------------------------------------------------
            base_day = obj.days(1);
            base_day_cells = find(obj.ds{base_day}.is_cell);
            base_day_num_cells = length(base_day_cells);
            
            % Scratch space for matching indices
            M = zeros(base_day_num_cells, num_days);
            M(:,1) = base_day_cells;
            
            % Fill out M columnwise
            for k = 2:num_days
                prev_day = obj.days(k-1);
                curr_day = obj.days(k);
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
        end % MultiDay
        
        % Accessors
        %------------------------------------------------------------
        function cell_idx = get_cell_idx(cell_idx, day_idx)
            cell_idx = obj.matched_indices(cell_idx, day_idx);
        end
    end
end
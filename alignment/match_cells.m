function matched_cells = match_cells(sources,m_base_to_others)
% Compute mutual matching matrix using pairwise cell overlap
% matrices
%
% inputs:
%   sources: 1 x [number of days] cell array containing paths to the
%       folders of the days of interest. The base day (the one to which all
%       other days are matched to) is taken to be the the first folder.
%   m_base_to_others: 1 x [number of days-1] cell array where each entry has
%       the pairwise overlap matrix between the base day and another.
%       m_base_to_others{k} must match with sources{k+1}.
%
% output:
%   matched_cells: [# of matched cells] x [# of days] array where each row 
%   contains the index of a single matched cell across the days of interest. 
%
% Hakan Inan (Mar 15)

num_days = length(sources);
classes = cell(1,num_days);

% Get class labels
for idx_day = 1:num_days
    class_day = load_classification(get_most_recent_file(sources{idx_day}, 'class_*.txt'));
    classes{idx_day} = strcmp(class_day,'not a cell') == 0;
end

num_ICs_base = length(classes{1}); % Use the 1st day as base
ICs_base = (1:num_ICs_base)';

% Match ICs
matched_ICs = zeros(num_ICs_base,num_days);
matched_ICs(:,1) = ICs_base;
idx_unmatched =[]; % Array of unmatched indices (w.r.t. base indices)
for idx_IC = ICs_base'
    is_IC_matched = 1;
    for idx_day = 2:num_days
        matched = m_base_to_others{idx_day-1}{idx_IC};
        if isempty(matched)
            matched_ICs(idx_IC,idx_day) = 0; 
            is_IC_matched = 0;
        else
            matched_ICs(idx_IC,idx_day) = matched(1); 
        end
    end
    
    if ~is_IC_matched % Call the IC unmatched 
        idx_unmatched = [idx_unmatched,idx_IC]; %#ok
    end    
end

matched_ICs(idx_unmatched,:) = []; % Remove unmatched rows
num_matched_ICs = size(matched_ICs,1);

% Update classes 
matched_classes = zeros(num_matched_ICs,num_days);
for idx_day = 1:num_days
    matched_classes(:,idx_day) = classes{idx_day}(matched_ICs(:,idx_day));
end

% Intersecion of the cell labels
idx_intrsct = find(prod(matched_classes,2)==1);

% Update cells
matched_cells = matched_ICs(idx_intrsct,:);

function match_list = load_matchlist(dataset_ids)
% Convenience function for making a 'match_list' to be used in the
% instantiation of a MultiDay object. The match results are loaded from
% pre-computed MAT files ('match_{X}_{Y}.mat') where {X} and {Y} are
% dataset IDs.

num_datasets = length(dataset_ids);
num_matches = num_datasets*(num_datasets-1)/2;
match_list = cell(num_matches, 4);

idx = 1;
for i = 1:num_datasets-1
    id_i = process_id(dataset_ids{i});
    for j = (i+1):num_datasets
        id_j = process_id(dataset_ids{j});
        match_data = load(sprintf('match_%s_%s.mat', id_i, id_j));
        
        match_list{idx,1} = id_i;
        match_list{idx,2} = id_j;
        match_list{idx,3} = match_data.m_1to2;
        match_list{idx,4} = match_data.m_2to1;
        
        idx = idx + 1;
    end
end

end % load_matchlist

function id = process_id(id)
    if isnumeric(id)
        id = num2str(id);
    end
end
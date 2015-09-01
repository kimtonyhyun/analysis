function [match_1to2, match_2to1, M] = match_masks(masks1, masks2)
% Computes the overlaps between two sets of logical masks. Output format is
% as follows:
%
% match_1to2{k} is a [Nx2] matrix where N is the number of masks from 
%   source 2 that match mask k of source 1. The first column of the matrix 
%   is the index of the matching mask in source 2; the second column is
%   second column the overlap measure between the two masks.
%

% Matching parameters
min_overlap_threshold = 1/3;

% Compute the matrix of mask overlaps
%------------------------------------------------------------
num_masks = [length(masks1) length(masks2)];
M = zeros(num_masks);
for i = 1:num_masks(1)
    if (mod(i,20)==0)
        fprintf('%s: Computing overlaps (%.1f%%)...\n',...
            datestr(now), 100*i/num_masks(1));
    end
    for j = 1:num_masks(2)
        M(i,j) = compute_mask_overlap(masks1{i}, masks2{j});
    end
end
fprintf('%s: Overlap matrix completed!\n', datestr(now));

% Find the nonzero elements of the mask overlap matrix
%------------------------------------------------------------
match_1to2 = find_nonzero_overlaps(M, min_overlap_threshold);
match_2to1 = find_nonzero_overlaps(M', min_overlap_threshold);  

end % match_masks

function match = find_nonzero_overlaps(M, threshold)

num_to_match = size(M,1);
match = cell(num_to_match,1);
for i = 1:num_to_match
    js = find(M(i,:) > threshold); % Col indices above threshold
    if isempty(js)
        match{i} = [];
    else
        overlaps = M(i,js);
        match{i} = [js' overlaps'];
        match{i} = sortrows(match{i}, 2); % Sort by overlap value
        match{i} = flipud(match{i}); % Descending
    end
end

end % find_nonzero_overlaps
function [match_1to2, match_2to1] = bijective_filter(match_1to2, match_2to1)
% Makes sure that the the matching is bijective, i.e.
%   1. Each cell in Dataset1 matches to at most one cell in Dataset 2 (and
%      visa versa)

for i = 1:length(match_1to2)
    match_itoj = match_1to2{i};
    num_matches = size(match_itoj,1);
    if (num_matches > 1)
        for k = 2:num_matches
            j = match_itoj(k,1);
            match_jtoi = match_2to1{j};
            [~, diff_inds] = setdiff(match_jtoi(:,1), i);
            if ~isempty(diff_inds)
                match_2to1{j} = match_jtoi(diff_inds,:);
            else
                match_2to1{j} = [];
            end
        end
        match_1to2{i} = match_itoj(1,:);
    end
end

for j = 1:length(match_2to1)
    match_jtoi = match_2to1{j};
    num_matches = size(match_jtoi,1);
    if (num_matches > 1)
        for k = 2:num_matches
            i = match_jtoi(k,1);
            match_itoj = match_1to2{i};
            [~, diff_inds] = setdiff(match_itoj(:,1), j);
            if ~isempty(diff_inds)
                match_1to2{i} = match_itoj(diff_inds,:);
            else
                match_1to2{i} = [];
            end
        end
        match_2to1{j} = match_jtoi(1,:);
    end
end

end % bijective_filter
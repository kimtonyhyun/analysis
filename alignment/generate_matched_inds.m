function matched_inds = generate_matched_inds(m_1to2)
% Convert a match output from 'run_alignment' between _two_ datasets into a
% list of matched indices. To generate matched indices involving more than
% two datasets, use the more sophisticated method via MultiDay instead.

num_max_matches = length(m_1to2);

matched_inds = zeros(num_max_matches, 2);

idx = 0;
for k = 1:num_max_matches
    if ~isempty(m_1to2{idx})
        idx = idx + 1;
        matched_inds(idx,:) = [idx m_1to2{idx}(1)];
    end
end
matched_inds = matched_inds(1:idx,:);

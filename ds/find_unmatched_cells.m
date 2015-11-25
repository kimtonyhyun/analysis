function [unmatched_ids, matched_ids] = find_unmatched_cells(ds1, m_1to2)
% Find classified cells in provided DaySummary ('ds1') with no match in
% the corresponding matching ('m_1to2')

is_cell = find(ds1.is_cell);

m_1to2 = m_1to2(is_cell);
is_not_matched = cellfun(@isempty, m_1to2);

unmatched_ids = is_cell(is_not_matched);
matched_ids = setdiff(is_cell, unmatched_ids);
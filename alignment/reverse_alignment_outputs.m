function [m_1to2_new, m_2to1_new, info_new] = reverse_alignment_outputs(m_1to2, m_2to1, info)
% Rearrange the outputs of 'run_alignment' as if the order of the
% DaySummary arguments were reversed.

m_1to2_new = m_2to1;
m_2to1_new = m_1to2;

info_new = info;

info_new.alignment.selected_cells = fliplr(info.alignment.selected_cells);
sc = info.alignment.selected_centers;
info_new.alignment.selected_centers = cat(3, sc(:,:,2), sc(:,:,1));
info_new.tform = invert(info.tform);
info_new.filters_2to1 = info.filters_1to2;
info_new.filters_1to2 = info.filters_2to1;
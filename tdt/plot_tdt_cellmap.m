function [A_pos, A_neg] = plot_tdt_cellmap(ds, tdt)
% Draw a cell map where tdTomato-positive neurons are shown in red, while
% tdTomato-negative neurons are shown in blue

A_pos = zeros(size(ds.cell_map_ref_img), 'single');
A_neg = A_pos;

num_pos_cells = 0;
for k = tdt.pos
    if ds.is_cell(k)
        num_pos_cells = num_pos_cells + 1;
        A_pos = A_pos + ds.cells(k).im;
    end
end

num_neg_cells = 0;
for k = tdt.neg
    if ds.is_cell(k)
        num_neg_cells = num_neg_cells + 1;
        A_neg = A_neg + ds.cells(k).im;
    end
end

C = A_pos - A_neg;

clims = prctile(abs(C(:)), 99) * [-1 1];
imagesc(C, clims);
axis image;
colormap redblue;
set(gca, 'TickLength', [0 0]);
title(sprintf('tdT-pos (%d cells; red) vs. -neg (%d cells; blue)',...
              num_pos_cells, num_neg_cells));
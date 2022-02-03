function [A_pos, A_neg] = plot_tdt_cellmap(ds, tdt)
% Draw a cell map where tdTomato-positive neurons are shown in red, while
% tdTomato-negative neurons are shown in blue

A_pos = zeros(size(ds.cell_map_ref_img), 'single');
A_neg = A_pos;

for k = tdt.pos
    A_pos = A_pos + ds.cells(k).im;
end

for k = tdt.neg
    A_neg = A_neg + ds.cells(k).im;
end

C = A_pos - A_neg;

clims = prctile(abs(C(:)), 99) * [-1 1];
imagesc(C, clims);
axis image;
colormap redblue;
set(gca, 'TickLength', [0 0]);
% title('tdTomato-positive (red) vs. -negative (blue)');
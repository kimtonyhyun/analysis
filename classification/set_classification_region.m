function set_classification_region(ds)
% Allows the user to graphically define a rectangular region in the FOV.
% All cell candidates outside of the region will be set to 'not a cell'.

hf = figure;
imagesc(ds.cell_map_ref_img);
axis image; colormap gray;
hold on;
for k = 1:ds.num_cells
    bd = ds.cells(k).boundary;
    plot(bd(:,1), bd(:,2), 'w');
end

fprintf('set_classification_roi: Please provide a rectangular region over the image.\n');
fprintf('  Double click on the rectangle when done.\n');
h_rect = imrect;
rect_params = round(wait(h_rect));
x_bounds = [rect_params(1) rect_params(1)+rect_params(3)-1];
y_bounds = [rect_params(2) rect_params(2)+rect_params(4)-1];
close(hf);

ds.reset_labels;
num_rejected = 0;
for k = 1:ds.num_cells
    com = ds.cells(k).com;
    if (com(1) < x_bounds(1)) || (com(1) > x_bounds(2)) ||...
       (com(2) < y_bounds(1)) || (com(2) > y_bounds(2))
       ds.cells(k).label = 'not a cell';
       num_rejected = num_rejected + 1;
    end
end

num_remaining = ds.num_cells - num_rejected;
cprintf('blue', 'set_classification_roi: %d cell candidates remaining\n', num_remaining);
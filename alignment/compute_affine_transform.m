function [info, masks1, masks2] = compute_affine_transform(ds1, ds2, alignment_cell_inds)
% IC map alignment based on user-defined control points.
%
% Inputs:
%   ds1/2: DaySummary object containing cell maps to be aligned
%   alignment_cell_inds: [N x 2] matrix containing the indices of matched
%       cells. Columns map to ds1 and ds2, respectively.
%
% Outputs:
%   info: Struct containing results of the affine transformation fit
%   masks1: Cell array containing thresholded masks for source 1
%   masks2: Cell array containing _transformed_ thresholded masks
%           for source 2. (NOTE: the image dimensions of masks2 match the
%           dimensions of source 1 masks!)
%

num_alignment_points = size(alignment_cell_inds, 1);

% Alignment is based on control _points_ (i.e. XY coordinates)
alignment_xy_coords = zeros(num_alignment_points, 2, 2); % XY position for each cell

% Historically, we used the mean of the boundary as the XY position for
% eahc cell. We don't use the filter's COM because the COM can deviate from
% the visual center of the cell's position as indicated by its boundary.
% The COM computed from the binarized mask may also be worth consideration.
for k = 1:num_alignment_points
    idx1 = alignment_cell_inds(k,1);
    bd1 = ds1.cells(idx1).boundary;
    alignment_xy_coords(k,:,1) = mean(bd1, 1);
    
    idx2 = alignment_cell_inds(k,2);
    bd2 = ds2.cells(idx2).boundary;
    alignment_xy_coords(k,:,2) = mean(bd2, 1);
end


% Transform Source2 onto Source1
%------------------------------------------------------------
tform = fitgeotrans(alignment_xy_coords(:,:,2),... % Moving points
                    alignment_xy_coords(:,:,1),... % Fixed points
                    'affine');

figure; % Pre-transform comparison
plot_boundaries_with_transform(ds1, 'b', 2, alignment_cell_inds(:,1));
hold on;
plot_boundaries_with_transform(ds2, 'r', 1, alignment_cell_inds(:,2));
hold off;
title('Pre-transform: Dataset1 (blue) vs. Dataset2 (red)');

figure; % Post-transform comparison
plot_boundaries_with_transform(ds1, 'b', 2, alignment_cell_inds(:,1));
hold on;
plot_boundaries_with_transform(ds2, 'r', 1, alignment_cell_inds(:,2), tform);
hold off;
title('Post-transform: Dataset1 (blue) vs. Dataset2 (red)');

% Prep output
%------------------------------------------------------------
masks1 = {ds1.cells.mask};
masks2 = {ds2.cells.mask};

mask1_ref = imref2d(size(masks1{1}));
for k = 1:ds2.num_cells
    masks2{k} = imwarp(masks2{k}, tform, 'OutputView', mask1_ref);
end

info.alignment.num_points = num_alignment_points;
info.alignment.selected_cells = alignment_cell_inds;
info.alignment.selected_centers = alignment_xy_coords;
info.tform = tform;

end % compute_affine_transform
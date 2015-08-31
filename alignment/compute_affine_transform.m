function [info, masks1, masks2] = compute_affine_transform(ds1, ds2)
% IC map alignment based on user-defined control points.
%
% Inputs:
%   ds1/2: DaySummary object containing cell maps to be aligned
%
% Outputs:
%   info: Struct containing results of the affine transformation fit
%   masks1: Cell array containing thresholded IC masks for source 1
%   masks2: Cell array containing _transformed_ thresholded IC masks
%           for source 2. (Note that the image dimensions match the
%           dimensions of source 1 masks!)
%

num_ics_for_alignment = 3;

% To programmatically address either of the two DaySummary's
ds = cell(1,2);
ds{1} = ds1;
ds{2} = ds2;

% Display the two sets of ICs
figure;
ax1 = subplot(121);
ds1.plot_cell_boundaries;
title('Dataset 1');

ax2 = subplot(122);
ds2.plot_cell_boundaries;
title('Dataset 2');

% Allow the user to select the ICs used in matching
%------------------------------------------------------------
sel_colors = 'ybm';
fprintf('compute_affine_transform: Please select %d cells from each dataset (in order)\n',...
    num_ics_for_alignment);

sel_ics = zeros(num_ics_for_alignment,2); % List of selected ICs
sel_ics_centers = zeros(num_ics_for_alignment, 2, 2); % XY position of each selected IC
sel_ics_idx = [0 0]; % Number of ICs selected from each dataset

% Loop until required number of ICs have been selected
while (~all(sel_ics_idx == num_ics_for_alignment))
    click_xy = round(ginput(1)); % Get user click
    if (gca == ax1) % Axis 1 was clicked
        source_idx = 1;
    elseif (gca == ax2)
        source_idx = 2;
    end
    
    ic_idx = ds{source_idx}.get_cell_by_xy(click_xy);
    if ~isempty(ic_idx) % Hit
        sel_idx = sel_ics_idx(source_idx) + 1;
        if (sel_idx <= num_ics_for_alignment)
            boundary = ds{source_idx}.cells(ic_idx).boundary;
            fill(boundary(:,1),...
                 boundary(:,2),...
                 sel_colors(mod(sel_idx,length(sel_colors))+1));

            sel_ic_centroid = mean(boundary,1);
            fprintf('  Dataset%d: Cell %d selected (at [%.1f %.1f])!\n',...
                source_idx, ic_idx,...
                sel_ic_centroid(1), sel_ic_centroid(2));
            
            sel_ics(sel_idx, source_idx) = ic_idx;
            sel_ics_idx(source_idx) = sel_idx;
            sel_ics_centers(sel_idx,:,source_idx) = sel_ic_centroid;
        else
            fprintf('  Dataset%d: No more ICs needed!\n',...
                source_idx);
        end
    else % No hit
        fprintf('  Dataset%d: No IC detected at cursor!\n',...
            source_idx);
    end
end
fprintf('  All reference ICs selected!\n');

% Transform Source2 onto Source1
%------------------------------------------------------------
tform = fitgeotrans(sel_ics_centers(:,:,2),... % Moving points
                    sel_ics_centers(:,:,1),... % Fixed points
                    'affine');

figure;
% subplot(121); % Pre-transform comparison
plot_boundaries(ds1, 'b', 2, sel_ics(:,1), []);
plot_boundaries(ds2, 'r', 1, sel_ics(:,2), []);
title('Pre-transform: Dataset1 (blue) vs. Dataset2 (red)');
axis equal;
set(gca, 'YDir', 'Reverse');

figure;
% subplot(122); % Post-transform comparison
plot_boundaries(ds1, 'b', 2, sel_ics(:,1), []);
plot_boundaries(ds2, c2, 'r', 1, sel_ics(:,2), tform);
title('Post-transform: Dataset1 (blue) vs. Dataset2 (red)');
axis equal;
set(gca, 'YDir', 'Reverse');

% Transform Source2 masks for output

mask1_ref = imref2d(size(masks{1}{1}));
for ic_idx = 1:ds2.num_cells
    masks{2}{ic_idx} = imwarp(masks{2}{ic_idx}, tform,...
        'OutputView', mask1_ref);
end

% Prep output
%------------------------------------------------------------
masks1 = masks{1};
masks2 = masks{2};

info.num_cells1 = ds1.num_cells;
info.num_cells2 = ds2.num_cells;
info.num_ics_for_alignment = num_ics_for_alignment;
info.ic_filter_threshold = ic_filter_threshold;
info.sel_ics = sel_ics;
info.sel_ics_centers  = sel_ics_centers;
info.tform = tform;

end % compute_affine_transform

function plot_boundaries(ds, linespec, linewidth, sel_ics, tform)
    % Plot boundaries as a single color, with an optional transform
    for k = 1:ds.num_cells
        boundary = ds.cells(k).boundary;
        if ~isempty(tform) % Optional spatial transform
            boundary = transformPointsForward(tform, boundary);
        end
        if ismember(k, sel_ics) % One of the user-selected ICs
            fill(boundary(:,1), boundary(:,2), linespec, 'LineWidth', linewidth);
        elseif ds.is_cell(k)
            plot(boundary(:,1), boundary(:,2), linespec, 'LineWidth', linewidth);
        end
        hold on;
    end
end

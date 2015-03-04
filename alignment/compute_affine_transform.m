function [info, masks1, masks2] = compute_affine_transform(source_dir1, source_dir2)
% IC map alignment based on user-defined control points.
%
% Inputs:
%   source_dirX: Path to directory containing ICA (ica_*.mat) and
%       classification (class_*.txt) outputs
%
% Outputs:
%   info: Struct containing results of the affine transformation fit
%   masks1: Cell array containing thresholded IC masks for source 1
%   masks2: Cell array containing _transformed_ thresholded IC masks
%           for source 2. (Note that the image dimensions match the
%           dimensions of source 1 masks!)
%

num_ics_for_alignment = 3;
ic_filter_threshold = 0.3;

% Load data
%------------------------------------------------------------
sources = {source_dir1, source_dir2};
s1 = load(get_most_recent_file(sources{1}, 'ica_*.mat'));
c1 = load_classification(get_most_recent_file(sources{1}, 'class_*.txt'));

s2 = load(get_most_recent_file(sources{2}, 'ica_*.mat'));
c2 = load_classification(get_most_recent_file(sources{2}, 'class_*.txt'));

bounds = cell(1,2); % Boundaries of ICs
masks  = cell(1,2); % BW image of thresholded ICs

% Display the two sets of ICs
figure;
ax1 = subplot(121);
[bounds{1}, masks{1}] = plot_ic_filters(s1.ica_filters, c1, ic_filter_threshold);
title(strrep(sources{1},'_','\_'));

ax2 = subplot(122);
[bounds{2}, masks{2}] = plot_ic_filters(s2.ica_filters, c2, ic_filter_threshold);
title(strrep(sources{2},'_','\_'));

% Allow the user to select the ICs used in matching
%------------------------------------------------------------
sel_colors = 'ybk';
fprintf('compute_affine_transform: Please select %d ICs from each dataset (in order)\n',...
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
    
    ic_idx = get_selected_ic_idx(bounds{source_idx}, click_xy);
    if ~isempty(ic_idx) % Hit
        sel_idx = sel_ics_idx(source_idx) + 1;
        if (sel_idx <= num_ics_for_alignment)
            boundary = bounds{source_idx}{ic_idx};
            fill(boundary(:,1),...
                 boundary(:,2),...
                 sel_colors(mod(sel_idx,length(sel_colors))+1));

            sel_ic_centroid = mean(boundary,1);
            fprintf('  %s: IC %d selected (at [%.1f %.1f])!\n',...
                sources{source_idx}, ic_idx,...
                sel_ic_centroid(1), sel_ic_centroid(2));
            
            sel_ics(sel_idx, source_idx) = ic_idx;
            sel_ics_idx(source_idx) = sel_idx;
            sel_ics_centers(sel_idx,:,source_idx) = sel_ic_centroid;
        else
            fprintf('  %s: No more ICs needed!\n',...
                sources{source_idx});
        end
    else % No hit
        fprintf('  %s: No IC detected at cursor!\n',...
            sources{source_idx});
    end
end
fprintf('  All reference ICs selected!\n');

% Transform Source2 onto Source1
%------------------------------------------------------------
tform = fitgeotrans(sel_ics_centers(:,:,2),... % Moving points
                    sel_ics_centers(:,:,1),... % Fixed points
                    'affine');

figure;
subplot(121); % Pre-transform comparison
plot_boundaries(bounds{1}, c1, 'b', 2, sel_ics(:,1), []);
plot_boundaries(bounds{2}, c2, 'r', 1, sel_ics(:,2), []);
title(sprintf('Pre-transform: %s (blue) vs. %s (red)',...
                strrep(sources{1}, '_', '\_'),...
                strrep(sources{2}, '_', '\_')));
axis equal;
set(gca, 'YDir', 'Reverse');

subplot(122); % Post-transform comparison
plot_boundaries(bounds{1}, c1, 'b', 2, sel_ics(:,1), []);
plot_boundaries(bounds{2}, c2, 'r', 1, sel_ics(:,2), tform);
title(sprintf('Post-transform: %s (blue) vs. %s (red)',...
                strrep(sources{1}, '_', '\_'),...
                strrep(sources{2}, '_', '\_')));
axis equal;
set(gca, 'YDir', 'Reverse');

% Transform Source2 masks for output
mask1_ref = imref2d(size(masks{1}{1}));
for ic_idx = 1:s2.ica_info.num_ICs
    masks{2}{ic_idx} = imwarp(masks{2}{ic_idx}, tform,...
        'OutputView', mask1_ref);
end

% Prep output
%------------------------------------------------------------
masks1 = masks{1};
masks2 = masks{2};

info.source1 = source_dir1;
info.source2 = source_dir2;
info.num_ics1 = s1.ica_info.num_ICs;
info.num_ics2 = s2.ica_info.num_ICs;
info.num_ics_for_alignment = num_ics_for_alignment;
info.ic_filter_threshold = ic_filter_threshold;
info.sel_ics = sel_ics;
info.sel_ics_centers  = sel_ics_centers;
info.tform = tform;

end % compute_affine_transform

function [boundaries, masks, ic_map] = plot_ic_filters(ica_filters, cl, ic_filter_threshold)
    % Plot grayscale image of IC filters along with their boundaries
    num_ics = size(ica_filters,3);
    boundaries = cell(num_ics, 1);
    masks = cell(num_ics, 1);
    
    ic_map = sum(ica_filters,3);
    imagesc(ic_map);
    colormap gray;
    axis equal;
    hold on;
    
    for ic_idx = 1:num_ics
        ic_filter = ica_filters(:,:,ic_idx);
        [ic_boundaries, mask] = compute_ic_boundary(ic_filter, ic_filter_threshold);
        ic_boundary = ic_boundaries{1}; % Longest boundary
        if any(strcmp(cl{ic_idx}, {'phase-sensitive cell', 'cell'}))
            color = 'g';
        else
            color = 'r';
        end
        plot(ic_boundary(:,1), ic_boundary(:,2), color);
        
        boundaries{ic_idx} = ic_boundary;
        masks{ic_idx} = mask;
    end
end

function plot_boundaries(boundaries, cl, linespec, linewidth, sel_ics, tform)
    % Plot boundaries as a single color, with an optional transform
    num_ics = length(cl);
    for ic_idx = 1:num_ics
        boundary = boundaries{ic_idx};
        if ~isempty(tform) % Optional spatial transform
            boundary = transformPointsForward(tform, boundary);
        end
        if ismember(ic_idx, sel_ics) % One of the user-selected ICs
            fill(boundary(:,1), boundary(:,2), linespec, 'LineWidth', linewidth);
        else
            if any(strcmp(cl{ic_idx}, {'phase-sensitive cell', 'cell'}))
                plot(boundary(:,1), boundary(:,2), linespec, 'LineWidth', linewidth);
            end
        end
        hold on;
    end
end

function selected_ic_idx = get_selected_ic_idx(boundaries, click_xy)
    % Return the first IC index whose filter boundary encloses the
    %   XY position provided in `click_xy`
    num_ics = length(boundaries);
    selected_ic_idx = []; % If no hit, then return empty
    for ic_idx = 1:num_ics
        bound = boundaries{ic_idx};
        if inpolygon(click_xy(1), click_xy(2), bound(:,1), bound(:,2))
            selected_ic_idx = ic_idx;
            break;
        end
    end
end
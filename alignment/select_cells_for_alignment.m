function selected_cell_inds = select_cells_for_alignment(ds1, ds2, num_alignment_points)

% To programmatically address either of the two DaySummary's
ds = cell(1,2);
ds{1} = ds1;
ds{2} = ds2;

% Display the two sets of ICs
figure;
ax1 = subplot(121);
ds1.plot_cell_boundaries('cells');
hold on;
title('Dataset 1');

ax2 = subplot(122);
ds2.plot_cell_boundaries('cells');
hold on;
title('Dataset 2');

% Allow the user to select the ICs used in matching
%------------------------------------------------------------
sel_colors = jet(num_alignment_points);
fprintf('Please select %d cells from each dataset (in order)\n', num_alignment_points);

selected_cell_inds = zeros(num_alignment_points, 2); % List of selected cells
num_selected = [0 0]; % Number of ICs selected from each dataset

% Loop until required number of ICs have been selected
while (~all(num_selected == num_alignment_points))
    click_xy = round(ginput(1)); % Get user click
    if (gca == ax1) % Axis 1 was clicked
        source_idx = 1;
    elseif (gca == ax2)
        source_idx = 2;
    end
    
    cell_idx = ds{source_idx}.get_cell_by_xy(click_xy, 'cells');
    if ~isempty(cell_idx) % Hit
        sel_idx = num_selected(source_idx) + 1;
        if (sel_idx <= num_alignment_points)
            boundary = ds{source_idx}.cells(cell_idx).boundary;
            fill(boundary(:,1), boundary(:,2), sel_colors(sel_idx,:));
            
            selected_cell_inds(sel_idx, source_idx) = cell_idx;
            num_selected(source_idx) = sel_idx;
            
            fprintf('  Dataset%d: Cell %d selected!\n', source_idx, cell_idx);
        else
            fprintf('  Dataset%d: No more cells needed!\n',...
                source_idx);
        end
    else % No hit
        fprintf('  Dataset%d: No cell detected at cursor!\n',...
            source_idx);
    end
end
fprintf('  All reference cells selected!\n');
function plot_boundaries(ds, varargin)
% Plot boundaries as a single color, with an optional transform. Can select
% a subset of cells to fill in.

line_color = 'b';
line_width = 1;
filled_cells = [];
tform = [];

% If 'show_all_cells' is false, we will show a limited number of boundaries
% sorted by their distance to 'center'. This is for performance reasons,
% as some DaySummary instances contain thousands of cells.
show_all_cells = true;
center = [0 0]';

for k = 1:length(varargin)
    if ischar(varargin{k})
        switch lower(varargin{k})
            case 'color'
                line_color = varargin{k+1};
            case {'linewidth', 'width'}
                line_width = varargin{k+1};
            case {'fill', 'filled_cells'}
                filled_cells = varargin{k+1};
            case {'transform', 'tform'}
                tform = varargin{k+1};
            case {'display_center', 'center'}
                show_all_cells = false;
                center = varargin{k+1};
                center = center(:); % Force column vector
        end
    end
end

cell_inds = find(ds.is_cell);
if ~show_all_cells
    num_to_show = min(50, ds.num_classified_cells);
    
    % We select subset to show based on distance to 'target' 
    coms = cell2mat({ds.cells(cell_inds).com}); % Each column is a COM
    if ~isempty(tform)
        coms = transformPointsForward(tform, coms')'; % Note: transformPointsForward operates on rows
    end
    distances = vecnorm(coms - center);
    [~, sorted_inds] = sort(distances, 'ascend');
    cell_inds = cell_inds(sorted_inds(1:num_to_show));
end

for k = cell_inds
    boundary = ds.cells(k).boundary;
    if ~isempty(tform) % Optional spatial transform
        boundary = transformPointsForward(tform, boundary);
    end
    if ismember(k, filled_cells)
        fill(boundary(:,1), boundary(:,2), line_color,...
             'LineWidth', line_width,...
             'FaceAlpha', 0.4);
    else
        plot(boundary(:,1), boundary(:,2), 'Color', line_color, 'LineWidth', line_width, 'HitTest', 'off');
    end
    hold on;
end
axis equal tight;
hold off;
set(gca, 'YDir', 'Reverse');
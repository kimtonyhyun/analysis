function plot_boundaries_with_transform(ds, linecolor, linewidth, filled_cells, tform, show_noncells)
    % Plot boundaries as a single color, with an optional transform. Can
    % subselect cells to be filled in.
    
    % TODO: Clean up this parameter parsing...
    if (nargin < 3)
        linewidth = 1;
        filled_cells = [];
        tform = [];
        show_noncells = false;
    elseif (nargin < 4)
        filled_cells = [];
        tform = [];
        show_noncells = false;
    elseif (nargin < 5)
        tform = [];
        show_noncells = false;
    elseif (nargin < 6)
        show_noncells = false;
    end   
    
    for k = 1:ds.num_cells
        boundary = ds.cells(k).boundary;
        if ~isempty(tform) % Optional spatial transform
            boundary = transformPointsForward(tform, boundary);
        end
        if ismember(k, filled_cells)
            fill(boundary(:,1), boundary(:,2), linecolor,...
                 'LineWidth', linewidth,...
                 'FaceAlpha', 0.4);
        else
            if ds.is_cell(k)
                plot(boundary(:,1), boundary(:,2), 'Color', linecolor, 'LineWidth', linewidth, 'HitTest', 'off');
            elseif show_noncells
                plot(boundary(:,1), boundary(:,2), ':', 'Color', linecolor, 'LineWidth', linewidth, 'HitTest', 'off');
            end
        end
        hold on;
    end
    axis equal tight;
    hold off;
    set(gca, 'YDir', 'Reverse');
end
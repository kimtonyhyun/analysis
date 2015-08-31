function plot_cell_boundaries(obj)
    % Display the cell map
    imagesc(obj.cell_map_ref_img);
    colormap gray;
    axis equal tight;

    hold on;
    for k = 1:obj.num_cells
        boundary = obj.cells(k).boundary;
        if obj.is_cell(k)
            color = 'g';
        else
            color = 'r';
        end
        plot(boundary(:,1), boundary(:,2), color);
    end
    hold off;
end
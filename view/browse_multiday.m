function browse_multiday(md)

h = figure;

common_cell_idx = 1;
while (1)
    draw_cell(common_cell_idx);
    
    % Ask user for command
    prompt = sprintf('Multiday browser (md cell %d of %d) >> ',...
        common_cell_idx, md.num_cells);
    resp = strtrim(input(prompt, 's'));
    
    val = str2double(resp);
    if (~isnan(val)) % Is a number
        
    else
        resp = lower(resp);
        switch (resp)
            case {'n', ''} % Next cell
                if (common_cell_idx < md.num_cells)
                    common_cell_idx = common_cell_idx + 1;
                else
                    fprintf('  Already at final md cell!\n');
                end
                
            case 'p' % Previous cell
                if (common_cell_idx > 1)
                    common_cell_idx = common_cell_idx - 1;
                else
                    fprintf('  Already at first md cell!\n');
                end
                
            case 'q' % Exit
                close(h);
                break;
                
            otherwise
                fprintf('  Could not parse "%s"\n', resp);
        end
    end
end

    function draw_cell(common_cell_idx)
        for k = 1:md.num_days
            day = md.valid_days(k);
            cell_idx_k = md.get_cell_idx(common_cell_idx, day);
            ds = md.day(day); % Just a shorthand
            
            % Draw image of cell on each day
            subplot(3, md.num_days, k);
            cell_image_k = ds.cells(cell_idx_k).im;
            cell_com_k = ds.cells(cell_idx_k).com;
            zoom_half_width = min(size(cell_image_k))/10;
            
            imagesc(cell_image_k);
            axis image;
            xlim(cell_com_k(1)+zoom_half_width*[-1 1]);
            ylim(cell_com_k(2)+zoom_half_width*[-1 1]);
            title(sprintf('Day %d -- Cell %d', day, cell_idx_k));
            
            % Draw raster
            subplot(3, md.num_days, [md.num_days+k 2*md.num_days+k]);
            ds.plot_cell_raster(cell_idx_k, 'draw_correct');
        end
    end % draw_cell

end % browse_multiday
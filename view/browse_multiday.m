function browse_multiday(md)

h = figure;

common_cell_idx = 1;
while (1)
    draw_md_cell(md, common_cell_idx);
    
    % Ask user for command
    prompt = sprintf('MD browser (Sort day %d, ID %d of %d) >> ',...
        md.sort_day, common_cell_idx, md.num_cells);
    resp = strtrim(input(prompt, 's'));
    
    val = str2double(resp);
    if (~isnan(val)) % Is a number
        if ((1 <= val) && (val <= md.num_cells))
            common_cell_idx = val;
        else
            fprintf('  Sorry, %d is not a valid md cell index\n', val);
        end
    else
        resp = lower(resp);
        if isempty(resp) % Empty string gets mapped to "n"
            resp = 'n';
        end
        
        switch resp(1)
            case 'd' % Day selected. View detailed raster on that day.
                % Parse the remainder of the string for day
                val = str2double(resp(2:end));
                if ~isnan(val)
                    if ismember(val, md.valid_days)
                        day_cell_idx = md.get_cell_idx(common_cell_idx, val);
                        view_detailed_raster(md.day(val), day_cell_idx);
                    else
                        fprintf('  Sorry, %d is not a valid day for this MultiDay\n', val);
                    end
                end
                
            case 's' % Sort the MD object on a different day.
                val = str2double(resp(2:end));
                if ~isnan(val)
                    if ismember(val, md.valid_days)
                        sort_inds = md.sort_matches_by_day(val);
                        % Keep the currently selected cell in view,
                        % despite the resorting of the underlying MD
                        common_cell_idx = find(sort_inds==common_cell_idx, 1);
                        fprintf('  MultiDay is now sorted on Day %d\n', val);
                    else
                        fprintf('  Sorry, %d is not a valid day for this MultiDay\n', val);
                    end
                end
                
            case 'c' % Jump to a cell on the current sorted day.
                val = str2double(resp(2:end));
                if ~isnan(val)
                    cells_from_sort_day = md.get_indices(md.sort_day);
                    common_cell_idx2 = find(cells_from_sort_day == val);
                    if ~isempty(common_cell_idx2)
                        common_cell_idx = common_cell_idx2;
                    else
                        fprintf('  Sorry, %d is not a valid cell on Day %d\n', val, md.sort_day);
                    end
                end
                
            case 't' % "Take" screenshot
                output_name = sprintf('md_%03d.png', common_cell_idx);
                print('-dpng', output_name);
                fprintf('  Screenshot saved to "%s"\n', output_name);
                
            case 'n'
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

end % browse_multiday
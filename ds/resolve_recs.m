function resolution_list = resolve_recs(md)

h = figure;

% For matched cell index k:
%   resolution_list(k,1): Selected rec index
%   resolution_list(k,2): Cell index within that rec
resolution_list = zeros(md.num_cells, 2);

common_cell_idx = 1;
while (1)
    draw_cell(common_cell_idx);
    
    % Ask user for command
    prompt = sprintf('Resolve recs (ID %d of %d) >> ',...
        common_cell_idx, md.num_cells);
    resp = strtrim(input(prompt, 's'));
    
    val = str2double(resp);
    if (~isnan(val)) % Is a number
        if ((1 <= val) && (val <= md.num_days)) % Check valid rec index
            if (md.matched_indices(common_cell_idx, val) ~= 0)
                resolution_list(common_cell_idx,1) = val;
                resolution_list(common_cell_idx,2) = ...
                    md.matched_indices(common_cell_idx, val);

                if (common_cell_idx < md.num_cells)
                    common_cell_idx = common_cell_idx + 1;
                else
                    fprintf('  Already at final md cell!\n');
                end
            else
                fprintf('  Sorry, %d is not a valid rec index for this cell\n', val);
            end
        else
            fprintf('  Sorry, %d is not a valid rec index\n', val);
        end
    else
        resp = lower(resp);
        if isempty(resp) % Empty string gets mapped to "n"
            resp = 'n';
        end
        
        switch resp(1)
            case 'f'
                common_cell_idx = 1;

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

    function draw_cell(common_cell_idx)
        clf;
        
        % Draw filters
        for k = 1:md.num_days
            day = md.valid_days(k);
            cell_idx_k = md.get_cell_idx(common_cell_idx, day);
            
            if (cell_idx_k ~= 0)
                % Draw filters
                ds_cell = md.day(day).cells(cell_idx_k);
                subplot(3, md.num_days, k);
                com = ds_cell.com;
                boundary = ds_cell.boundary;
                
                imagesc(ds_cell.im);
                colormap gray;
                axis image;
                hold on;
                plot(boundary(:,1),boundary(:,2),'Color',get_color(k),...
                     'LineWidth', 2);
                hold off;
                zoom_half_width = min(size(ds_cell.im))/20;
                xlim(com(1)+zoom_half_width*[-1 1]);
                ylim(com(2)+zoom_half_width*[-1 1]);
                
                title_str = sprintf('Iter %d -- Cell %d', day, cell_idx_k);
                title(title_str);
            end
        end
        
        % Draw traces       
        h_trace = subplot(3, md.num_days, (md.num_days+1):(3*md.num_days));
        set(zoom(h_trace), 'Motion', 'horizontal');
        trace_offset = 0;
        for k = fliplr(1:md.num_days)
            day = md.valid_days(k);
            cell_idx_k = md.get_cell_idx(common_cell_idx, day);
            
            if (cell_idx_k ~= 0)
                trace_k = md.day(day).get_trace(cell_idx_k);
                plot(trace_k + trace_offset, 'Color', get_color(k));
                trace_offset = trace_offset + max(trace_k);
                hold on;
            end
        end
        xlim([0 length(trace_k)]);
        ylim([0 trace_offset]);
        grid on;
        xlabel('Frame');
        ylabel('Traces');
        set(gca, 'YTickLabel', []);

        function color = get_color(day_idx)
            colors = {'b', 'r', [0 0.5 0]};
            color = colors{mod(day_idx-1, md.num_days)+1};
        end 
        
    end % draw_cell

end % resolve_recs
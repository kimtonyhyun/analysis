function res_list = resolve_merged_recs(md, M, varargin)
% Used to "resolve" duplicate cell filters for a given dataset. Used, for
% example, during iterative cell sorting or generating "union" filters
% (i.e. transferring unmatched filters across sessions).
%
% Input 'md' is a MultiDay instance created by 'create_merge_md'.
%
% For matched cell index k, 'res_list' indicates:
%   res_list(k,1): Selected rec index.
%                  Special flags: 0 (unassigned), -1 (omit)
%   res_list(k,2): Cell index within that rec

res_list = zeros(md.num_cells, 2);

% Defaults
if isempty(M)
    state.clim = [0 1];
else
    state.clim = compute_movie_scale(M);
end
normalize_traces = false;
rec_names = cell(md.num_days, 1);
for i = 1:md.num_days
    rec_names{i} = sprintf('Rec %d', i);
end
rec_colors = flipud(winter(md.num_days));

for i = 1:length(varargin)
    vararg = varargin{i};
    if ischar(vararg)
        switch lower(vararg)
            case {'norm_traces', 'normalize_traces'}
                normalize_traces = true;
            case 'names'
                rec_names = cellfun(...
                    @(x) strrep(x, '_', '\_'),... % Escape underscores in names
                    varargin{i+1},...
                    'UniformOutput', false);
                assert(length(rec_names) == md.num_days,...
                    'Provided list of rec names does not match number of datasets in MultiDay');
        end
    end
end

% Set up GUI
%------------------------------------------------------------
h_fig = figure;
set(h_fig, 'DefaultAxesTitleFontWeight', 'normal');
set(h_fig, 'WindowButtonUpFcn', @end_drag);
set(h_fig, 'WindowScrollWheelFcn', @scroll_frame);

h_recs = cell(md.num_days, 1);
h_rec_ims = cell(md.num_days, 1);
for i = 1:md.num_days
    h_recs{i} = subplot(3, md.num_days, i);
    colormap gray;
end

h_trace_sp = subplot(3, md.num_days, (md.num_days+1):(3*md.num_days));
h_traces = cell(md.num_days, 1);
hold on;
num_frames = md.day(1).full_num_frames;
dummy_trace_vals = -Inf*ones(1,num_frames);
for i = 1:md.num_days
    h_traces{i} = plot(1:num_frames, dummy_trace_vals,...
        'Color', rec_colors(i,:),...
        'HitTest', 'off');
end
set(zoom(h_trace_sp), 'Motion', 'horizontal');
set(h_trace_sp, 'YTick', []);
set(h_trace_sp, 'TickLength', [0 0]);
xlabel('Frame');
ylabel('Traces');
xlim([1 num_frames]);
h_time = plot(-1*[1 1], [0 Inf], 'k', 'ButtonDownFcn', @start_drag); % Time indicator
hold off;
set(h_trace_sp, 'ButtonDownFcn', @click_handler);

% Auxiliary display parameters
[height, width] = size(md.day(1).cell_map_ref_img);
zoom_half_width = min([height width])/20;
trace_zoom_half_width = num_frames/20;

% For convenience, if a filter only shows up on one rec,
% then assign it automatically
%------------------------------------------------------------
for j = 1:md.num_cells
    inds = md.matched_indices(j,:);
    num_recs = sum(inds>0);
    if (num_recs == 1)
        rec_ind = find(inds, 1);
        res_list(j,:) = [rec_ind inds(rec_ind)];
    end
end


cell_idx = 1;
last_drawn_idx = 0;

while (1)
    if (last_drawn_idx ~= cell_idx)
        draw_match(cell_idx);
        state.zoomed = false;
        last_drawn_idx = cell_idx;
    end
    
    % Ask user for command
    %------------------------------------------------------------
    switch res_list(cell_idx,1)
        case 0
            assign_str = 'unassigned';
        case -1
            assign_str = 'omit';
        otherwise
            assign_str = rec_names{res_list(cell_idx,1)};
    end
    prompt = sprintf('Resolve recs (ID %d of %d: %s) >> ',...
        cell_idx, md.num_cells, assign_str);
    resp = strtrim(input(prompt, 's'));
    
    val = str2double(resp);
    if (~isnan(val)) % Is a number
        if ((1 <= val) && (val <= md.num_days)) % Check valid rec index
            if (md.matched_indices(cell_idx, val) ~= 0)
                res_list(cell_idx,1) = val;
                res_list(cell_idx,2) = ...
                    md.matched_indices(cell_idx, val);

                go_to_next_cell();
            else
                fprintf('  Sorry, %d is not a valid rec index for this cell\n', val);
            end
        else
            fprintf('  Sorry, %d is not a valid rec index\n', val);
        end
    else
        resp = lower(resp);
        if isempty(resp)
            go_to_next_cell();
        else
            switch resp(1)
                case 'c' % Jump to a cell
                    val = str2double(resp(2:end));
                    if ~isnan(val)
                        if ((1 <= val) && (val <= md.num_cells))
                            cell_idx = val;
                        else
                            fprintf('  Sorry, %d is not a valid cell on Day %d\n', val, md.sort_day);
                        end
                    end

                case 'f' % "filter"
                    % Force redraw of the GUI, which re-displays the filter
                    last_drawn_idx = 0;

                case {'h', 'l'} % "higher/lower contrast"
                    c_range = diff(state.clim);
                    if (strcmp(resp, 'h'))
                        new_clim = state.clim + c_range*[0.1 -0.1];
                        fprintf('  Increased contrast (new CLim=[%.3f %.3f])\n',...
                            new_clim(1), new_clim(2));
                    else
                        new_clim = state.clim + c_range*[-0.1 0.1];
                        fprintf('  Decreased contrast (new CLim=[%.3f %.3f])\n',...
                            new_clim(1), new_clim(2));
                    end

                    for i = 1:md.num_days
                        h_sp = h_recs{i};
                        if ~isempty(h_sp)
                            set(h_sp, 'CLim', new_clim);
                        end
                    end
                    state.clim = new_clim;

                case {'n', 'x'} % Mark cell to be omitted
                    res_list(cell_idx,1) = -1;
                    go_to_next_cell();

                case 'p' % Previous cell
                    if (cell_idx > 1)
                        cell_idx = cell_idx - 1;
                    else
                        fprintf('  Already at first md cell!\n');
                    end

                case 'z'
                    if state.zoomed
                        xlim(h_trace_sp, [1 num_frames]);
                    else
                        % Get the current frame from the time indicator
                        cf = get(h_time, 'XData'); cf = cf(1);
                        xlim(h_trace_sp, cf + trace_zoom_half_width * [-1 +1]);
                    end
                    state.zoomed = ~state.zoomed;

                case 'q' % Exit
                    close(h_fig);
                    break;

                otherwise
                    fprintf('  Could not parse "%s"\n', resp);
            end
        end
    end 
end

    function go_to_next_cell()
        to_be_assigned = (res_list(:,1) == 0);
        next_idx = find_next_cell_to_process(cell_idx, to_be_assigned);
        if isempty(next_idx)
            cprintf([0 0.5 0], '  All cells have been assigned!\n');
        else
            cell_idx = next_idx;
        end
    end

    function draw_match(common_cell_idx)        
        title(h_trace_sp, '');
        set(h_time, 'XData', [0 0]); % Hide the cursor
        
        trace_offset = 0;
        for k = fliplr(1:md.num_days) % For stacking of traces
            day = md.valid_days(k);
            cell_idx_k = md.get_cell_idx(common_cell_idx, day);

            if (cell_idx_k == 0) % No matching cell in this Rec
                blank_subplot(h_recs{k});
                h_rec_ims{k} = [];
                set(h_traces{k}, 'YData', dummy_trace_vals);
            else
                % Draw filters
                %-------------------------------------------------------
                ds_cell = md.day(day).cells(cell_idx_k);
                com = ds_cell.com;
                boundary = ds_cell.boundary;

                subplot(h_recs{k});
                h_rec_ims{k} = imagesc(rescale_filter_to_clim(ds_cell.im, state.clim), state.clim);
                hold on;
                plot(boundary(:,1), boundary(:,2),...
                     'Color', rec_colors(k,:),...
                     'LineWidth', 2);
                hold off;
                axis image;
                xlim(com(1)+zoom_half_width*[-1 1]);
                ylim(com(2)+zoom_half_width*[-1 1]);
                
                title_str = sprintf('\\bf(%d)\n\\rm%s, cell %d',...
                    day, rec_names{k}, cell_idx_k);
                title(title_str);
                
                % Draw traces
                %-------------------------------------------------------
                trace_k = md.day(day).get_trace(cell_idx_k);
                if normalize_traces
                    trace_k_min = min(trace_k);
                    trace_k_max = max(trace_k);
                    trace_k = (trace_k - trace_k_min)/(trace_k_max - trace_k_min);
                end
                set(h_traces{k}, 'YData', trace_k  + trace_offset);
                trace_offset = trace_offset + max(trace_k);
            end
        end
        xlim(h_trace_sp, [1 num_frames]);
        ylim(h_trace_sp, [0 trace_offset]);
        set(h_time, 'YData', [0 trace_offset]);
        
    end % draw_cell

    function blank_subplot(h_sp)
        subplot(h_sp);
        plot([1 width], [height 1], 'k');
        hold on;
        plot([1 width], [1 height], 'k');
        hold off;
        title('');
        xlim([1 width]);
        ylim([1 height]);
        set(h_sp, 'XTick', []);
        set(h_sp, 'YTick', []);
        axis image;
    end % blank_subplot

    % Functions required for interactive display
    %-------------------------------------------------------
    function click_handler(~, e)
        k = round(e.IntersectionPoint(1));
        render_frame(k);
    end

    function start_drag(~, ~)
        set(h_fig, 'WindowButtonMotionFcn', @drag_frame);
    end

    function end_drag(~, ~)
        set(h_fig, 'WindowButtonMotionFcn', '');
    end

    function drag_frame(~, ~)
        cp = get(h_trace_sp, 'CurrentPoint');
        k = round(cp(1));
        render_frame(k);
    end

    function scroll_frame(~, e)
        k = get(h_time, 'XData'); k = k(1);
        if (e.VerticalScrollCount < 0) % Scroll up
            k = k - 1;
        else
            k = k + 1;
        end
        render_frame(k);
    end

    function render_frame(k)
        k = max(1,k); k = min(k, num_frames); % Clamp
        set(h_time, 'XData', k*[1 1]);
        title(h_trace_sp, sprintf('Frame %d', k));
        
        if ~isempty(M)
            A = M(:,:,k);
            
            for k = 1:md.num_days
                if ~isempty(h_rec_ims{k})
                    set(h_rec_ims{k}, 'CData', A);
                end
            end
        end
    end % render_frame)

end % resolve_recs
function [raster, trial_inds] = plot_cell_raster(obj, cell_idx, varargin)
% Plots a raster of cell activity, where trials are aligned to the closing
% gate frame.
%
% Optional argument 'draw_correct' will place a box at the end
% of each trial indicating correct (green) or incorrect (red)
%
% Additional optional arguments allow for filtering of trials,
% e.g. "plot_cell_raster(cell_idx, 'start', 'east')"
%

    align_index = 3; % By default, align to the closing of gate
    
    display_trial = ones(obj.num_trials, 1);
    draw_correct = 0;
    
    if ~isempty(varargin)
        for k = 1:length(varargin)
            vararg = varargin{k};
            if ischar(vararg)
                switch lower(vararg)
                    case 'align'
                        align_index = varargin{k+1};
                        
                    case 'draw_correct'
                        draw_correct = 1;
                end
            end
        end
        % Trial filtering arguments
        display_trial = obj.filter_trials(varargin{:});
    end

    num_filtered_trials = sum(display_trial);

    switch align_index
        case 1
            align_str = 'Frames relative to trial start';
        case 2
            align_str = 'Frames relative to gate open';
        case 3
            align_str = 'Frames relative to gate close';
        case 4
            align_str = 'Frames relative to trial end';
    end
    [pre_offset, post_offset] = compute_frame_offsets(obj.trial_indices, align_index);
    num_trunc_frames = post_offset-pre_offset+1;

    raster = zeros(num_filtered_trials, num_trunc_frames);
    correctness = zeros(num_filtered_trials);

    counter = 0;
    for k = 1:obj.num_trials
        if display_trial(k)
            counter = counter + 1;

            % Compute aligned frame indices into current trial
            ti = obj.trial_indices(k,:);
            ti = ti - (ti(1)-1);
            af = ti(align_index); % alignment frame
            pre_frame = af + pre_offset;
            post_frame = af + post_offset;

            tr = obj.get_trace(cell_idx, k);
            raster(counter,:) = tr(pre_frame:post_frame);

            correctness(counter) = obj.trials(k).correct;
        end
    end

    imagesc(pre_offset:post_offset, 1:num_filtered_trials, raster, 'HitTest', 'off');
    colormap jet;
    xlabel(align_str);
    ylabel('Trial index');
    set(gca, 'TickLength', [0 0]);

    if (draw_correct)
        corr_width = 0.025*num_trunc_frames;
        for k = 1:num_filtered_trials
            if correctness(k)
                corr_color = 'g';
            else
                corr_color = 'r';
            end
            rectangle('Position', [post_offset k-0.5 corr_width 1],...
                      'FaceColor', corr_color);
        end
        xlim([pre_offset post_offset+corr_width]);
    end
    
    trial_inds = find(display_trial);
end
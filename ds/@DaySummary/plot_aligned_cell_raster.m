function raster = plot_aligned_cell_raster(obj, cell_idx, varargin)
% Optional argument 'draw_correct' will place a box at the end
% of each trial indicating correct (green) or incorrect (red)
%
% Additional optional arguments allow for filtering of trials,
% e.g. "plot_cell_raster(cell_idx, 'start', 'east')"
%
% Note: May merge into 'plot_cell_raster' in the future...

    display_trial = ones(obj.num_trials, 1);
    draw_correct = 0;
    if ~isempty(varargin)
        for k = 1:length(varargin)
            vararg = varargin{k};
            if ischar(vararg)
                switch lower(vararg)
                    case 'draw_correct'
                        draw_correct = 1;
                end
            end
        end
        % Trial filtering arguments
        display_trial = obj.filter_trials(varargin{:});
    end

    num_filtered_trials = sum(display_trial);

    [pre_offset, post_offset] = compute_offsets(obj.trial_indices);
    num_trunc_frames = post_offset-pre_offset+1;

    raster = zeros(num_filtered_trials, num_trunc_frames);
    correctness = zeros(num_filtered_trials);

    counter = 0;
    for k = 1:obj.num_trials
        if display_trial(k)
            counter = counter + 1;

            ti = obj.trial_indices(k,:);
            ti = ti - (ti(1)-1);
            cgf = ti(3);
            pre_frame = cgf + pre_offset;
            post_frame = cgf + post_offset;

            tr = obj.get_trace(cell_idx, k);
            raster(counter,:) = tr(pre_frame:post_frame);

            correctness(counter) = obj.trials(k).correct;
        end
    end

    imagesc(pre_offset:post_offset, 1:num_filtered_trials, raster);
    colormap jet;
    xlabel('Frames relative to gate close');
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
end

function [pre_offset, post_offset] = compute_offsets(frame_indices)
    close_gate_frames = frame_indices(:,3);
    frame_indices = frame_indices - repmat(close_gate_frames, [1 4]); % Align each trial to closing of gate
    
    pre_offset = max(frame_indices(:,1));
    post_offset = min(frame_indices(:,4));
end
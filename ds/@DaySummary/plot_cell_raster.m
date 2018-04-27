function [raster, info] = plot_cell_raster(obj, cell_idx, varargin)
% Plots a raster of cell activity, where trials are aligned to the closing
% gate frame by default. Other alignment points are also possible via the 
% 'align' vararg. Additional arguments allow for filtering of trials, e.g.
%
%   plot_cell_raster(cell_idx, 'start', 'east')
%
% See 'DaySummary.get_aligned_trace' for full details.
%
% Optional argument 'draw_correct' will place a box at the end
% of each trial indicating correct (green) or incorrect (red)
%
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
    end
   
    % Raster: [num_cells x num_aligned_frames]
    [raster, info] = obj.get_aligned_trace(cell_idx, varargin{:});
    
    switch info.align_idx
        case 1
            align_str = 'Frames relative to trial start';
        case 2
            align_str = 'Frames relative to gate open';
        case 3
            align_str = 'Frames relative to gate close';
        case 4
            align_str = 'Frames relative to trial end';
    end
    
    imagesc(info.aligned_time, 1:info.num_trials, raster, 'HitTest', 'off');
    colormap jet;
    xlabel(align_str);
    ylabel('Trial index');
    set(gca, 'TickLength', [0 0]);

    if (draw_correct)
        corr_width = 0.025*size(raster,2);
        trial_inds = info.trial_inds;
        for k = 1:length(trial_inds)
            trial_idx = trial_inds(k);
            if obj.trials(trial_idx).correct
                corr_color = 'g';
            else
                corr_color = 'r';
            end
            rectangle('Position', [info.aligned_time(end) k-0.5 corr_width 1],...
                      'FaceColor', corr_color);
        end
        xlim([info.aligned_time(1) info.aligned_time(end)+corr_width]);
    end
end
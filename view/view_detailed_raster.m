function view_detailed_raster(ds, cell_idx)
% Plot a detailed "raster summary" of a single cell in DaySummary
%
% Usage:
%     m3d15 = DaySummary(sources.maze, 'rec001', 'reconst', 'excludeprobe');
%       14-Aug-2015 13:58:23: Loaded data from rec001/rec_150813-115003.mat
%       14-Aug-2015 13:58:23: Loaded classification from rec001/class_150813-115320.txt
%     browse_rasters(m3d15);
%

while (1)
    clf;
    draw_rasters(); 
    
    % Ask user for command
    prompt = sprintf('Cell %d raster >> ', cell_idx);
    resp = strtrim(input(prompt, 's'));
    
    val = str2double(resp);
    if (~isnan(val)) % Is a number (trial index)
        if ~ds.is_behavior_loaded
            fprintf('  Behavior video not loaded into DaySummary!\n');
        else
            if ((1 <= val) && (val <= ds.num_trials))
                fprintf('  Showing trial %d. Press any key to return.\n', val);
                draw_trial(val);
            else
                fprintf('  Error, %d is not a valid trial index!\n', val);
            end
        end
    else
        resp = lower(resp);
        switch (resp)
            case 't' % "Take" screenshot
                savename = sprintf('raster_%d.png', cell_idx);
                print('-dpng', savename);
                fprintf('  Saved raster as "%s"!\n', savename);
            
            case {'q', ''} % Exit
                break;
                
            otherwise
                fprintf('  Could not parse "%s"\n', resp);
        end
    end
end

    % Helper functions
    %------------------------------------------------------------
    function draw_rasters()
        % Image of cell
        subplot(3,4,[1 2]);
        imagesc(ds.cells(cell_idx).im);
        axis image;
        title(sprintf('Cell %d (%s)', cell_idx, ds.cells(cell_idx).label));
        colormap jet; freezeColors;
        
        % Raster of all trials, with correctness
        subplot(3,4,[5 6 9 10]);
        ds.plot_cell_raster(cell_idx, 'draw_correct');
        raster_scale = get(gca, 'CLim'); % Scale that applies to all trials
        title('All trials');
        colormap jet; freezeColors;
        
        % Divide rasters by correctness
        subplot(3,4,3);
        ds.plot_cell_raster(cell_idx, 'correct');
        set(gca, 'CLim', raster_scale);
        title('Correct');
        subplot(3,4,4);
        ds.plot_cell_raster(cell_idx, 'incorrect');
        set(gca, 'CLim', raster_scale);
        title('Incorrect');
        
        % Divide rasters by start location
        subplot(3,4,7);
        ds.plot_cell_raster(cell_idx, 'start', 'west');
        set(gca, 'CLim', raster_scale);
        title('West start'); 
        subplot(3,4,8);
        ds.plot_cell_raster(cell_idx, 'start', 'east');
        set(gca, 'CLim', raster_scale);
        title('East start');
        
        % Divide rasters by end location
        subplot(3,4,11);
        ds.plot_cell_raster(cell_idx, 'end', 'south');
        set(gca, 'CLim', raster_scale);
        title('South end');
        subplot(3,4,12);
        ds.plot_cell_raster(cell_idx, 'end', 'north');
        set(gca, 'CLim', raster_scale);
        title('North end'); 
    end % draw_rasters

    function draw_trial(trial_idx)
        Mb = ds.get_behavior_trial(trial_idx); % Behavior movie
        tr = ds.trials(trial_idx).traces(cell_idx,:);
        num_frames_in_trial = length(tr);
        
        % Show trace
        subplot(3,4,[3 4]);
        plot(1:num_frames_in_trial, tr);
        xlim([1 num_frames_in_trial]);
        y_min = min(tr);
        y_max = max(tr);
        y_range = y_max - y_min;
        ylim([y_min y_max] + 0.1*y_range*[-1 1]);
        grid on;
        title(sprintf('Trial %d', trial_idx));
        xlabel('Frames');
        ylabel('Signal [a.u.]');
        
        % Show behavior movie
        subplot(3,4,[7 8 11 12]);
        imagesc(Mb(:,:,1));
        set(gca, 'XTick', []);
        set(gca, 'YTick', []);
        colormap gray;
        pause;
    end % draw_trial
end % view_cell_rasters
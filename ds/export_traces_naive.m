function [X, xs, ys] = export_traces_naive(md, trial_map, extent)
% X = EXPORT_TRACES(md, trial_map)
%
% Exports traces of all trials specified by trial_map into a 3D matrix 'X'
%   with dimensions [neurons x time x trials]. Each trial has been time
%   normalized to have the same length (in number of samples).
%
% Also provides the mouse's trajectory on each trial. Trajectories are also
%   time normalized.
%

    num_cells = md.num_cells;
    num_samples = compute_resample_size(md); % Num samples per trial
    num_trials = size(trial_map,1);
    
    X = zeros(num_cells, num_samples, num_trials);
    
    xs = cell(num_trials,1);
    ys = cell(num_trials,1);
    
    resample_grid = linspace(0, 1, num_samples);
    for k = 1:num_trials
        % day and neuron indices
        di = trial_map(k,1);
        ni = md.get_indices(di);

        % traces for this trial
        trial = md.day(di).trials(trial_map(k,2));
        [t,x,y] = truncate_trial(trial.centroids, trial.start, extent);
        traces = trial.traces(ni,t);
        
        % Resample each trial to common number of samples
        orig_grid = linspace(0, 1, size(traces,2));
        X(:,:,k) = interp1(orig_grid, traces', resample_grid)';
        xs{k} = interp1(orig_grid, x, resample_grid);
        ys{k} = interp1(orig_grid, y, resample_grid);
    end

    function [t_idx,x,y] = truncate_trial(xy, start_arm, extent)
        % x and y coordinates
        x = xy(:,1);
        y = xy(:,2);

        % truncate as specified
        if strcmp(extent, 'full')
            t_idx = true(size(x));   
        elseif strcmp(extent, 'first')
            % NOTE (THK): I'm concerned about the 'extent' handling, since
            % it is possible for the mouse to leave and come back into the
            % "pass zone" which means that we could be selecting a
            % discontinuous chunk of a trial. This sort of selection could
            % have weird behavior when temporally warping. Thus, it would
            % be preferrable to make it so that only a single chunk is
            % retrieved by 'extent' handling.
            %
            % Additionally, the hard coded parameters may not work for all
            % datasets. We took data from two separate experimental setups
            % and the position of the camera wasn't identical.
            switch start_arm
                case 'east'
                    t_idx = y > east_start_boundary(x);
                case 'west'
                    t_idx = y < west_start_boundary(x);
                otherwise
                    error('extent not implemented for probe trials');
            end
        elseif strcmp(extent,'second')
            switch start_arm
                case 'east'
                    t_idx = y < east_start_boundary(x);
                case 'west'
                    t_idx = y > west_start_boundary(x);
                otherwise
                    error('extent not implemented for probe trials');
            end
        else
            error('extent not specified correctly')
        end

        x = x(t_idx,:);
        y = y(t_idx,:);

    end % truncate_trial

    function y = east_start_boundary(x)
        y = -x+600;
    end

    function y = west_start_boundary(x)
        y = -x+450;
    end
end % export_traces_naive

function num_samples = compute_resample_size(md)
% Compute the common number of samples to use when resampling traces.
    trial_indices = [];
    for di = md.valid_days
        trial_indices = [trial_indices; md.day(di).trial_indices]; %#ok<AGROW>
    end
    num_frames_per_trial = trial_indices(:,4) - trial_indices(:,1) + 1;

    num_samples = max(num_frames_per_trial);
end % compute_resample_size
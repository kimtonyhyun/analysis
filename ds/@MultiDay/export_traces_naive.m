function [X,x,y] = export_traces_naive(md, trial_map, extent)
% X = EXPORT_TRACES(md, trial_map)
%
% Exports traces of all trials specified by trial_map into
% a cell array X.

    num_trials = size(trial_map,1);
    X = cell(num_trials,1);
    x = cell(num_trials,1);
    y = cell(num_trials,1);
    for k = 1:num_trials
        % day and neuron indices
        d = trial_map(k,1);
        ni = md.matched_indices(:,md.valid_days == d);

        % traces for this trial
        trial = md.day(d).trials(trial_map(k,2));
        [t,x{k},y{k}] = truncate_trial(trial.centroids, trial.start, extent);
        X{k} = trial.traces(ni,t);
    end

    function [t_idx,x,y] = truncate_trial(xy, start_arm, extent)
        % x and y coordinates
        x = xy(:,1);
        y = xy(:,2);

        % truncate as specified
        if strcmp(extent, 'full')
            t_idx = true(size(x));   
        elseif strcmp(extent, 'first')
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
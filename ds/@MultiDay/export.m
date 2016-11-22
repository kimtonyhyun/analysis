function [X, meta, neuron_map, trial_map] = export(md, varargin)
% [X, meta, neuron_map, trial_map] = EXPORT(md)
%
% Exports cross-day aligned cell traces (via MultiDay) into a lightweight
% format (description below) for tensor analysis
%
% Output format:
%   X: Column cell vector [M x 1] where M is the total number 
%      of trials across all days of the requested 'trial_type' (e.g. 'en')
%       
%      X{j} is a matrix of cell traces, in the format
%      [num_matched_cells x num_frames_in_trial_j].
%
%   neuron_map: Matrix [num_matched_cells x num_days] that maps each 
%      matched neuron to its per-day neuron index.
%
%   trial_map: Matrix [M x 2] where M is the total number of trials of 
%      requested 'trial_type'. For the i-th trial, K(i,1) is the day index 
%      of that trial, and K(i,2) is the trial index in the original day.

    % By default, export full traces. Allow user to specify
    % otherwise with the 'extent' optional argument.
    extent = 'full';
    for k = 1:length(varargin)-1
        v = varargin{k};
        if ischar(v) && strcmp('extent',v)
            extent = varargin{k+1};
        end
    end

    % get trial map (filtering out those specified)
    trial_map = filter_trials(md, varargin{:});

    % neuron map for matching cells across days
    neuron_map = md.matched_indices;

    % activity traces for each trial
    [X,x,y] = export_traces(md, trial_map, extent);

    % metadata for each trial (start, end, turn, correct, etc.)
    meta = export_metadata(md, trial_map);
    meta.x = x; % also export position
    meta.y = y;

end % export

function meta = export_metadata(md, trial_map)
% X = EXPORT_METADATA(md, trial_map)
%
% Exports traces of all trials specified by trial_map into
% a cell array X.

    num_trials = size(trial_map,1);

    % get moving average of turn probability on each day
    ndays = length(md.valid_days);
    tp = cell(ndays);
    for di = 1:ndays
        d = md.valid_days(di);
        tp{di} = est_turn_probabilities(md.day(d));
    end

    % copy selected trials into lightweight cell array
    meta.start = cell(num_trials,1);
    meta.end = cell(num_trials,1);
    meta.correct = zeros(num_trials,1);
    meta.day = zeros(num_trials,1);
    meta.turn = cell(num_trials,1);
    meta.turn_prob = zeros(num_trials,1);
    for k = 1:num_trials
        % day and neuron indices
        d = trial_map(k,1);
        trial = md.day(d).trials(trial_map(k,2));

        % basic metadata associated with this trial
        meta.start{k} = trial.start;
        meta.end{k} = trial.end;
        meta.correct(k) = trial.correct;
        meta.day(k) = d;
        meta.turn{k} = trial.turn;

        % turn probability estimated for this trial
        di = md.valid_days == d;
        meta.turn_prob(k) = tp{di}(trial_map(k,2));
    end

    % mark each trial as allo vs ego-centric
    meta.strategy = cell(num_trials,1);
    e0 = NaN; % trial index of last east start
    w0 = NaN; % trial index of last west start
    for k = 1:num_trials
        pk = meta.turn_prob(k);
        if strcmp(meta.start{k},'east')
            if ~isnan(w0) && ~isnan(pk)
                pl = meta.turn_prob(w0);
                if pk > 0.99 && pl > 0.99
                    meta.strategy{k} = 'ego-right';
                elseif pk < 0.01 && pl < 0.01
                    meta.strategy{k} = 'ego-left';
                elseif pk > 0.99 && pl < 0.01
                    meta.strategy{k} = 'allo-north';
                elseif pk < 0.01 && pl > 0.99
                    meta.strategy{k} = 'allo-south';
                end
            end
            e0 = k;
        elseif strcmp(meta.start{k},'west')
            if ~isnan(e0) && ~isnan(pk)
                pl = meta.turn_prob(e0);
                if pk > 0.99 && pl > 0.99
                    meta.strategy{k} = 'ego-right';
                elseif pk < 0.01 && pl < 0.01
                    meta.strategy{k} = 'ego-left';
                elseif pk > 0.99 && pl < 0.01
                    meta.strategy{k} = 'allo-south';
                elseif pk < 0.01 && pl > 0.99
                    meta.strategy{k} = 'allo-north';
                end
            end
            w0 = k;
        else
            meta.strategy{k} = 'probe';
        end
        if isempty(meta.strategy{k})
            meta.strategy{k} = 'NA';
        end
    end
end % export_metadata

function [X,x,y] = export_traces(md, trial_map, extent)
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
end % export_traces
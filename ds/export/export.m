function [X, meta, neuron_map, trial_map] = export(md, varargin)
% [X, meta, neuron_map, trial_map] = EXPORT(md)
%
% Exports cross-day aligned cell traces from MultiDay into a lightweight
% format for further analysis (e.g. tensor analysis, decoding, etc.)
%
% Output format:
%
%   X: Cross-day traces formatted into a [neurons x time x trials] matrix.
%       Note that all trials will be formatted to the same length.
%
%   meta: Metadata associated with each trial (e.g. start location, etc.)
%
%   neuron_map: Matrix [neurons x num_days] that maps each 
%      matched neuron to its per-day neuron index.
%
%   trial_map: Matrix [trials x 2] where trial_map(i,1) is the day index 
%      and trial_map(i,2) is the trial index in the original day.
%

    timewarp_method = 'naive';

    extent = 'full'; % Used with timewarp_method == 'naive'
    align_idx = 3; % Used with timewarp_method == 'align'
    
    for k = 1:length(varargin)
        vararg = varargin{k};
        if ischar(vararg)
            switch lower(vararg)
                case 'extent'
                    extent = varargin{k+1}; % 'full', 'first', or 'second'
            end
        end
    end
    
    % get trial map (filtering out those specified)
    trial_map = md.filter_trials(varargin{:});

    % neuron map for matching cells across days
    neuron_map = md.matched_indices;

    % activity traces for each trial
    switch timewarp_method
        case 'naive'
            % 'naive' method will take each trial, and rescale time so
            % that each trial has the same number of samples.
            fprintf('  Exporting MD with NAIVE time warping method (extent="%s")...\n', extent);
            [X,x,y] = export_traces_naive(md, trial_map, extent);
            
        case 'align'
            % 'align' method will align each trial to one of four
            % intra-trial events:
            %   (1) Start of trial
            %   (2) Opening of gate
            %   (3) Closing of gate
            %   (4) End of trial
            % taking as many samples that is common to all trials (i.e.
            % bounded by the shortest trial). Time is not scaled, and can
            % be directly compared between trials.
            fprintf('  Exporting MD with ALIGN time warping method (align_idx=%d)...\n', align_idx);
            [X,x,y] = export_traces_align(md, trial_map, align_idx);
            
        otherwise
            error('Time warping method not recognized.');
    end
    
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
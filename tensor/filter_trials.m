function trial_map = filter_trials(dataobj,varargin)
% trial_map = FILTER_TRIALS(md,opts...)
%
%     Given a MultiDay or DaySummary object, return a trial_map specifying
%     the trials that satisfied the specified criteria. The following
%     options are available:
%       *  'start'   : any subset of {'north','south','east','west'}
%       *  'end'     : any subset of {'north','south','east','west'}
%       *  'correct' : any subset of {0,1}

    if isa(dataobj,'MultiDay')
        trial_map = filter_multiday(dataobj,varargin{:});
    elseif isa(dataobj,'DaySummary')
        trial_map = filter_day(dataobj,varargin{:});
    else
        error('input to filter_trials needs to be MultiDay or DaySummary object')
    end

end

function trial_map = filter_multiday(md,varargin)
    
    % preallocate
    max_trials = 0;
    for d = md.valid_days
        max_trials = max_trials + length(md.day(d).trials);
    end
    trial_map = zeros(max_trials,2);
    
    % filter each day
    i = 1;
    for d = md.valid_days
        day_map = filter_day(md.day(d),varargin{:});
        j = i-1 + length(day_map);
        trial_map(i:j,1) = d;
        trial_map(i:j,2) = day_map;
        i = j+1;
    end
    
    trial_map = trial_map(1:j,:);
end

function trial_idx = filter_day(ds,varargin)

    % parse optional inputs
    p = inputParser;
    p.addParameter('start', {'east','west'});
    p.addParameter('end', {'north','south','east','west'});
    p.addParameter('correct', [0,1]);
    p.parse(varargin{:});
    res = p.Results;

    K = length(ds.trials);
    trial_idx = zeros(K,1);
    i = 0;
    for k = 1:K
        trial = ds.trials(k);
        s = any(strcmp(trial.start,res.start));
        e = any(strcmp(trial.end,res.end));
        c = any(trial.correct == res.correct);
        if s && e && c
            i = i+1;
            trial_idx(i) = k;
        end
    end
    trial_idx = trial_idx(1:i,:);
    
end

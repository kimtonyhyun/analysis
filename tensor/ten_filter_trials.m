function [i1,i2,labels] = ten_filter_trials(session,trial_idx,type)
% [i1,i2,labels] = FILTER_TRIALS(session,trial_idx)
%
%     Returns indices ie and iw so that trial_idx(ie) contains east trials
%     and trial_idx(iw) contains west trials 

% synonyms
if strcmp(type,'error')
    type = 'correct';
elseif strcmp(type,'stop')
    type = 'end';
end

switch type
    case 'start'
        s1 = 'east';
        s2 = 'west';
        labels = {s1, s2};
    case 'end'
        s1 = 'north';
        s2 = 'south';
        labels = {s1, s2};
    case 'correct'
        s1 = '1';
        s2 = '0';
        labels = {'correct','error'};
    otherwise
        error(['unsupported trial filter ',type])
end

% east and west starts
%#ok<*AGROW>
i1 = [];
i2 = [];
for it = 1:length(trial_idx)
    trial = session.trials(trial_idx(it));
    trialdata = trial.(type);
    if isnumeric(trialdata)
        trialdata = num2str(trialdata);
    end
    if strcmp(trialdata,s1)
        i1 = [i1; it];
        continue
    end
    if strcmp(trialdata,s2)
        i2 = [i2; it];
        continue
    end
    error('assertion error')
end

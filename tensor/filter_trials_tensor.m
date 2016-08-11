function filtered_trial_indices = filter_trials_tensor(ds, trial_type)
% Return the trial indices in DaySummary that correspond to a trial type:
%
% Trial types:
%   'all': All non-probe trials
%   'east': All east starts
%   'west': All west starts
%   'en': All east -> north trials
%   'es': All east -> south trials
%   'wn': All west -> north trials
%   'ws': All west -> south trials
%   'everything': All trials (including probe)
%
% Note: This function is largely redundant with the built-in method
% 'filter_trials' of DaySummary. Leaving it for now, and will clean up
% later (THK).

all_trials = 1:ds.num_trials;
[en, es, wn, ws, probes] = unique_paths(ds);

switch trial_type
    case 'east'
        filtered_trial_indices = union(en,es);
    case 'west'
        filtered_trial_indices = union(wn,ws);
    case 'en'
        filtered_trial_indices = en;
    case 'es'
        filtered_trial_indices = es;
    case 'wn'
        filtered_trial_indices = wn;
    case 'ws'
        filtered_trial_indices = ws;
    case 'all' % remove probes
        filtered_trial_indices = setdiff(all_trials, probes);
    case 'everything' % keep probes
        filtered_trial_indices = all_trials;
    otherwise
        error('trial type not understood')
end

filtered_trial_indices = filtered_trial_indices'; % Return as a ROW vector
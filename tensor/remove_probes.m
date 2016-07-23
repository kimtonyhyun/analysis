function [idx] = remove_probes(session)
% [idx] = REMOVE_PROBES(session)
%
%     Given a session summary, return trial indices of non-probe trials.

idx = [];

for it = 1:session.num_trials
	s = session.trials(it).start;
	if strcmp('east',s) || strcmp('west',s)
		idx = [idx;it];
	end
end

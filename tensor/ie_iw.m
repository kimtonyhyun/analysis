function [ie,iw] = ie_iw(session,trial_idx)

% east and west starts
ie = [];
iw = [];
for it = 1:length(trial_idx)
	trial = session.trials(trial_idx(it));
	if strcmp(trial.start,'east')
		ie = [ie;it];
		continue
	end
	if strcmp(trial.start,'west')
		iw = [iw;it];
		continue
	end
	error('assertion error')
end

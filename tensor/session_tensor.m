function [X,cell_idx,trial_idx] = session_tensor(session)
% SESSION_TENSOR, construct a 3d tensor (neurons x time x trial) from a behavioral session.
%
%     X,cell_idx = SESSION_TENSOR(session)

% find cells
cell_idx = [];
for a = 1:session.num_cells
	if strcmp('not a cell',session.cells(a).label)
		continue
	end
	cell_idx = [cell_idx; a];
end

% determine median trial duration
lens = zeros(session.num_trials,1);
for b = 1:session.num_trials
	trial = session.trials(b);
	lens(b) = size(trial.traces,2);
end
len = ceil(median(lens));

% remove probe trials
trial_idx = remove_probes(session);

% linear interpolation to align trials.
X = zeros(length(cell_idx),len,length(trial_idx));
for b = 1:length(trial_idx)
	it = trial_idx(b);
	trial = session.trials(it);
	L = size(trial.traces,2);
	d = 1;
    for c = transpose(cell_idx)
        X(d,:,b) = interp1(linspace(1,len,L),trial.traces(c,:),1:len);
        d = d+1;
    end
end

% TO DO: Time warping? Position mapping?

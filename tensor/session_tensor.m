function [X,idx] = session_tensor(ds)
% SESSION_TENSOR, construct a 3d tensor (neurons x time x trial) from a behavioral session.
%
%     X,idx = SESSION_TENSOR(day_summary)

% find cells
idx = [];
for a = 1:ds.num_cells
	if strcmp('not a cell',ds.cells(a).label)
		continue
	end
	idx = [idx; a];
end

lens = zeros(ds.num_trials,1);
for b = 1:ds.num_trials
	trial = ds.trials(b);
	lens(b) = size(trial.traces,2);
end
len = ceil(median(lens));

X = zeros(length(idx),len,ds.num_trials);
for b = 1:ds.num_trials
	trial = ds.trials(b);
	L = size(trial.traces,2);
	d = 1;
    for c = transpose(idx)
        X(d,:,b) = interp1(linspace(1,len,L),trial.traces(c,:),1:len);
        d = d+1;
    end
end

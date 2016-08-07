function [X,cell_idx,trial_idx] = session_tensor(session,trial_type)
% SESSION_TENSOR, construct a 3d tensor (neurons x time x trial) from a behavioral session.
%
%     X,cell_idx = SESSION_TENSOR(session)
%     X,cell_idx = SESSION_TENSOR(session,'all') % all non-probe trials (default)
%     X,cell_idx = SESSION_TENSOR(session,'east') % all east starts
%     X,cell_idx = SESSION_TENSOR(session,'west') % all west starts
%     X,cell_idx = SESSION_TENSOR(session,'en') % all east -> north trials
%     X,cell_idx = SESSION_TENSOR(session,'es') % all east -> south trials
%     X,cell_idx = SESSION_TENSOR(session,'wn') % all west -> north trials
%     X,cell_idx = SESSION_TENSOR(session,'ws') % all west -> south trials

if nargin == 1
    trial_type = 'all';
end

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

% select trials of interest
everything = 1:length(session.trials);
[en,es,wn,ws,probes] = unique_paths(session);
switch trial_type
    case 'east'
        trial_idx = union(en,es);
        disp('building tensor to analyze east trials')
    case 'west'
        trial_idx = union(wn,ws);
        disp('building tensor to analyze west trials')
    case 'en'
        trial_idx = en;
        disp('building tensor to analyze east -> north trials')
    case 'es'
        trial_idx = es;
        disp('building tensor to analyze east -> south trials')
    case 'wn'
        trial_idx = wn;
        disp('building tensor to analyze west -> north trials')
    case 'ws'
        trial_idx = ws;
        disp('building tensor to analyze west -> south trials')
    case 'all' % remove probes
        trial_idx = setdiff(everything,probes);
        disp('building tensor to analyze all non-probe trials')
    case 'everything' % keep probes
        trial_idx = everything;
        disp('building tensor to analyze everything (including probes)')
    otherwise
        error('trial type not understood')
end

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

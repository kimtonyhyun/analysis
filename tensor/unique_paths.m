function [en,es,wn,ws,probes] = unique_paths(session)
% [en,es,wn,ws,probes] = UNIQUE_PATHS(session)
%
%     Given a session summary, return four lists of trial indices, which
%     indicate which of the four unique paths the mouse took on each trial.
%       *  `en` contains the "east -> north" trials
%       *  `es` contains the "east -> south" trials
%       *  `wn` contains the "west -> north" trials
%       *  `ws` contains the "west -> south" trials
%       *  `probes` contains all trials starting in north/south

%#ok<*AGROW>
en = [];
es = [];
wn = [];
ws = [];
probes = [];

for it = 1:session.num_trials
	x = session.trials(it).start;
    y = session.trials(it).end;
    if strcmp('north',x) || strcmp('south',x)
        probes = [probes;it];
    elseif strcmp('east',x)
        if strcmp('north',y)
            en = [en;it];
        elseif strcmp('south',y)
            es = [es;it];
        else
            error('invalid trial')
        end
    elseif strcmp('west',x)
        if strcmp('north',y)
            wn = [wn;it];
        elseif strcmp('south',y)
            ws = [ws;it];
        else
            error('invalid trial')
        end
    else
        error('invalid trial')
    end
end

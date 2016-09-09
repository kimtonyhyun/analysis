function trial_colors = get_trial_colors(md,trial_map,trial_type)
% trial_colors = GET_TRIAL_COLORS(md, trial_map, trial_type)
%
% Returns a K x 3 matrix with each row a rgb color mapped based on the
% trial_type parameter:
%
% GET_TRIAL_COLORS(md, trial_map, 'start')   % color based on start location
% GET_TRIAL_COLORS(md, trial_map, 'end')     % color based on end location
% GET_TRIAL_COLORS(md, trial_map, 'correct') % color based on correctness

% trial coloring labels
red = [1.0 0.0 0.0];
blue = [0.0 0.7 1.0];
switch trial_type
    case 'start'
        tc = {'east', 'E', blue;
              'west', 'W', red};
    case 'end'
        tc = {'north', 'N', blue;
              'south', 'S', red};
    case 'correct'
        tc = {'1', '1', blue;
              '0', '0', red};
    case 'strategy'
        error('This still needs to be implemented')
    otherwise
        error('Unknown trial_type parameter');
end

% set trial colors
nk = size(trial_map,1);
trial_colors = zeros(nk,3);
if ~isempty(tc)
    for k = 1:nk
        trial = md.day(trial_map(k,1)).trials(trial_map(k,2));
        trialdata = num2str(trial.(trial_type));
        idx = strcmp(tc(:,1), trialdata);
        trial_colors(k,:) = tc{idx,3};
    end
end

function [ num_persev_errors, num_random_errors ] = determine_error_type( file, strategy_switch, trial_start, trial_end)
% Reads behavior txt file and given pre/post switch strategies,
%       returns # of perseverative errors
% inputs:
%   strategy switch: elan, elas,eran,eras,anel,asel, aner, aser
%   trial_start, trial_end: first_tria, last_trial, e.g. 56, 155
%   file: plus maze text file
%
% parse plusmaze
[~,maze_data,~] = parse_plusmaze(file);

% consider only trials specified
maze_data = maze_data(trial_start:trial_end,:);

    start_arm = maze_data(:,1);
    end_arm = maze_data(:,2);
    actual_end_arm = maze_data(:,3);
    

% for given strategy switch, designate switching arm (where persev errors
% occur)
if strcmp(strategy_switch,'elan')
    persev_start = 'east';
elseif strcmp(strategy_switch,'anel')
    persev_start = 'east';
elseif strcmp(strategy_switch,'eras')
    persev_start = 'east';
elseif strcmp(strategy_switch,'aser')
    persev_start = 'east';
elseif strcmp(strategy_switch,'elas')
    persev_start = 'west';
elseif strcmp(strategy_switch,'asel')
    persev_start = 'west';
elseif strcmp(strategy_switch,'eran')
    persev_start = 'west';
elseif strcmp(strategy_switch,'aner')
    persev_start = 'west';
end

% find incorrect trials
incorr_trials = find(~strcmp(end_arm, actual_end_arm));

% count how many incorrect trials have a start in perseverative path
num_persev_errors = sum(strcmp(start_arm(incorr_trials),persev_start));

% assign all other errors as random errors
num_random_errors = size(incorr_trials,1)-num_persev_errors;

end


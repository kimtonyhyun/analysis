function [ percent_correct ] = calc_behavior_score( file , varargin )
% Reads behavior text file and gives behavior score as % correct
% file: expecting plusmaze txt file
% varargin:
% * 'trials': specify certain trials to consider, e.g.
%           e.g.: perc = % calc_behavior_score(file, 'trials', [1 55])
%                will return perc correct for trials 1 through 55 only
%           * expecting a matrix of trials to consider
% Updated 2015-03-12 Fori Wang

    trial_chunk = 0;
    if ~isempty(varargin)
            num_vararg = length(varargin);
            for k = 1:num_vararg
                switch varargin{k}
                    case 'trials' %specify trials to consider
                        trial_chunk = 1;
%                         trials_to_consider = varargin{k+1};
                         trial_start = varargin{k+1};
                         trial_end = varargin{k+2};
                        fprintf('Calculating score for trials given. \n');
%                      case 'skip' %specify trials (e.g. probe) not to consider
%                          skip_trials = 1;
%                          trials_to_skip = vararg{k+1}; 
                end
            end
    end

    % grab maze data from plusmaze file
    [~,maze_data,~] = parse_plusmaze(file);
    num_trials = size(maze_data,1);
    
    % if trials aren't specified, consider entire session
    if ~trial_chunk
%         trials_to_consider = 1:num_trials;
        trial_start = 1;
        trial_end = num_trials;
    end
    
    % determine correct choice and actual choice made by mouse
    end_arm = maze_data(trial_start:trial_end,2);
    actual_end_arm = maze_data(trial_start:trial_end,3);
    
    % calculate perfect correct
    chunk_length = size(end_arm,1);
    num_correct_instances = sum(strcmp(end_arm,actual_end_arm));
    percent_correct = num_correct_instances*100/chunk_length;

end
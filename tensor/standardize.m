function [dest] = standardize(X,trial_map)
% STANDARDIZE(X, trial_map), standardizes a set of multiday
% traces 
%
%     Args
%     ----
%     X = STANDARDIZE(X,trial_map,[method='soft'])

% default standardization
if nargin == 2
    method = 'soft';
end

% dimensions
n_cell = size(X{1},1);
n_trial = length(X);
n_time = ones(1,n_trial+1);
for k = 1:n_trial
    n_time(k+1) = size(X{k},2);
end

% preallocate final dest
dest = cell(size(X));

% list of days
days = unique(trial_map(:,1));

% standardization routines
if strcmp(method,'soft')
    % matricize full dataset
    x = zeros(n_cell,sum(n_time)-1);
    for k = 1:n_trial
        a = n_time(k);
        b = n_time(k+1)+a-1;
        x(:,a:b) = X{k};
    end
    % standardize each neuron
    denom = range(x,2);
    for c = 1:size(X,1)
        for k = 1:n_trial
            dest{k}(c,:) = X{k}(c,:) ./ denom(c);
        end
    end
elseif strcmp(method,'medmax')
    % max fluorescence for each cell on each trial
    mx = zeros(n_cell,n_trial);
    for k = 1:n_trial
        mx(:,k) = max(X{k},[],2);
    end
    % standardize within day
    for d = days' 
        % median max fluorescence on day d
        kd = find(trial_map(:,1)==d);
        med_X = median(mx(:,kd),2); % consider max here.
        MAX = max(med_X);
        
        % normalize each cell on day d
        for k = kd'
            dest{k} = zeros(size(X{k}));
            for i = 1:n_cell
                dest{k}(i,:) = X{k}(i,:) ./ ( 2*med_X(i) - MAX);
            end
        end
    end
else
    error('Standardization method not recognized.')
end

function [cmx] = cellmax(X, trial_map)

% separate days/sessions for this dataset
days = unique(trial_map(:,1));

K = length(X);     % num trials
N = size(X{1},1);  % num neurons
D = length(days);  % num days

cmx = zeros(N,D);

for k = 1:K
    % day for this trial
    d = find(days == trial_map(k,1));

    % max for each cell on this trial
    mx = max(X{k},[],2);

    % reduce op across all trials
    cmx(:,d) = max([ cmx(:,d) mx], [], 2);
end

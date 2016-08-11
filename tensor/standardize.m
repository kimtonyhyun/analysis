function [dest] = standardize(X,trial_map,method)
% SOFT_NORMALIZE, apply soft-normalization to each neuron in data tensor X and
% return a new, normalized tensor Y.
%
%     X = STANDARDIZE(X,trial_map,[method='soft'])

if nargin == 2
    method = 'soft';
end

n_cell = size(X{1},1);
n_trial = length(X);

dest = cell(size(X));
day = trial_map(:,1);
trial = trial_map(:,2);
days = unique(trial_map(:,1));

if strcmp(method,'soft')
    % max fluorescence for each cell on each trial
    mx = zeros(n_cell,n_trial);
    for k = 1:n_trial
        mx(:,k) = max(X{k},[],2);
    end
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
function [dest] = timewarp(X,method)
% TIMEWARP, given a cell array of trials, warp trials to achieve a
% consistent time axis. The result is a 3d tensor (neurons x time x
% trials).
%
%     X = TIMEWARP(X,[method='naive'])

if nargin == 1
    method = 'naive'; % default warping method
end

n_cell = size(X{1},1);
n_trial = length(X);

if strcmp(method,'naive')
    %% linear interpolation to align trials. %%
    
    % determine length
    lens = zeros(n_trial,1);
    for i = 1:n_trial
        lens(i) = size(X{i},2);
    end
    n_time = ceil(median(lens));
    
    % do interpolation.
    dest = zeros(n_cell,n_time,n_trial);
    for k = 1:n_trial
        old_ax = linspace(1,n_time,size(X{k},2));
        new_ax = 1:n_time;
        for i = 1:n_cell
            dest(i,:,k) = interp1(old_ax,X{k}(i,:),new_ax);
        end
    end
else
    error('Warping method not recognized.')
end

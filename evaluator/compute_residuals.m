function [Resids, Resids_n] = compute_residuals(X, Xest)

assert(all(size(X)==size(Xest)),...
    'Raw and fitted data do not have the same dimensions!');

Resids = squeeze(sum((X-Xest).^2,2)); % [num_neurons x num_trials]
Resids_n = sum(Resids,2);
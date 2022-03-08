function d = compute_trial_raster_dimensionality(trial_raster, method)
% Given a trial raster [num_trials x time_within_trial] for a single cell,
% compute the "dimensionality" of the cell's activity.

switch (method)
    case {'dim', 'dimensionality'}
        % Williams et al. Neuron (2018), Eq. 16
        C = trial_raster * trial_raster'; % [num_trials x num_trials]
        evs = eig(C); % Eigenvalues
        d = sum(evs)^2 / sum(evs.^2);
        
    case 'erank'
        % Roy and Vetterli. EURASIP 2007
        % Note: We normalize the result using the maximum possible rank
        [M, N] = size(trial_raster);
        Q = min([M, N]); % Maximum possible rank

        sigs = svd(trial_raster); % Singular values in descending order
        S = sum(sigs);
        p = sigs/S;
        H = -sum(p.*log(p));
        d = exp(H)/Q;
        
    otherwise
        error('Unrecognized method "%s"', method);
end

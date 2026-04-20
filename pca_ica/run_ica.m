function run_ica(pca_source, num_ICs, mu)
% Usage:
%   run_ica('_data/pca_n500.mat', 200, 0.1); % pca_source, num_ICs, mu

% Load PCA results
pca = load(pca_source);

% ICA
%------------------------------------------------------------
fprintf('%s: Computing ICA weights...\n', datestr(now));

ica_mixed = compute_spatiotemporal_ica_input(pca.filters, pca.traces, mu);

term_tol = 1e-5; % Termination tolerance
max_iter = 750;  % Max iterations of FastICA
ica_W = compute_ica_weights(ica_mixed, num_ICs, term_tol, max_iter)';

info.type = 'ica';
info.pca_source = pca_source;
info.num_pairs = num_ICs;
info.mu = mu; %#ok<*STRNU>

% Save ICA results
%------------------------------------------------------------
timestamp = datestr(now, 'yymmdd-HHMMSS');
icaw_savename = sprintf('icaw_%s.mat', timestamp);
save(icaw_savename, 'info', 'ica_W'); % Save weights

fprintf('%s: Computing ICA pairs (filters, traces) from weights...\n', datestr(now));
[filters, traces] = compute_ica_pairs(pca_source, icaw_savename);
save_rec(info, filters, traces);

fprintf('%s: All done!\n', datestr(now));
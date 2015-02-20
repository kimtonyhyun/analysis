function run_ica(pca_source, num_ICs, mu)

% Usage:
%   run_ica('_data/pca_n500.mat', 200, 0.1); % pca_source, num_ICs, mu

% Load PCA results
load(pca_source, 'pca_info', 'pca_filters', 'pca_traces');

% ICA
%------------------------------------------------------------
fprintf('%s: Computing ICA weights...\n', datestr(now));

ica_mixed = compute_spatiotemporal_ica_input(pca_filters, pca_traces, mu);

term_tol = 1e-5; % Termination tolerance
max_iter = 750;  % Max iterations of FastICA
ica_W = compute_ica_weights(ica_mixed, num_ICs, term_tol, max_iter)'; %#ok<NASGU>

ica_info.pca_source = pca_source;
ica_info.num_ICs = num_ICs;
ica_info.mu = mu; %#ok<*STRNU>

% Save ICA results
%------------------------------------------------------------
timestamp = datestr(now, 'yymmdd-HHMMSS');
ica_savename  = sprintf('ica_%s.mat',  timestamp);
icaw_savename = sprintf('icaw_%s.mat', timestamp);

save(icaw_savename, 'ica_info', 'ica_W'); % Save weights

fprintf('%s: Computing ICA pairs (filters, traces) from weights...\n', datestr(now));
[ica_filters, ica_traces] = compute_ica_pairs(pca_source, icaw_savename); %#ok<NASGU,ASGLU>
save(ica_savename,  'ica_info', 'ica_traces', 'ica_filters'); % Save ICA pairs

fprintf('%s: All done!\n', datestr(now));
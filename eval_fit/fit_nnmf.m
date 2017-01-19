function [N, T, Xest] = fit_nnmf(X, r)
% Fit NNMF model of rank r by unwrapping data tensor X.
%
% Output format:
%   N(:,i) = Neuron loadings of the i-th factor
%   T(:,j,i) = Temporal loading of the i-th factor, j-th trial
%   Xest: Reconstructed tensor, same dimension as the data tensor X
%

[num_neurons, num_samples, num_trials] = size(X);
Xu = reshape(X, [num_neurons, num_samples*num_trials]); % Unwrap along time

options = statset('nnmf');
options.Display = 'iter'; % Verbose output
options.MaxIter = 200;

% 'nnmf' output format:
%   N(:,i) = Neuron loadings of the i-th factor
%   T(i,:) = Temporal (unwrapped) loading of the i-th factor
[N,T] = nnmf(Xu, r, 'options', options,...
             'algorithm', 'als');
                
Xest = N*T; % Reconstructed tensor
Xest = reshape(Xest, [num_neurons, num_samples, num_trials]); % Rewrap into tensor

T = reshape(T, [r, num_samples, num_trials]);
T = permute(T, [2 3 1]); % [Intra-trial time x trials x rank]
function [A, B, Xest] = fit_nnmf(X, unfold_dim, r)
% Fit NNMF model of rank r by unfolding the data tensor X along a specified
% dimension.
%
% Output format:
%   - A(:,i) = Loadings on the unfolded dimension of the i-th factor
%   - B(:,j,i) =
%       For neural unfolding: Samples of the i-th factor, j-th TRIAL
%       For trial unfolding: Samples of the i-th factor, j-th NEURON
%   - Xest: Reconstructed tensor, same dimension as the data tensor X
%
% Usage examples:
%   - Unfold on neurons: [N, TrialTime, Xest] = fit_nnmf(X, 1, r);
%   - Unfold on trials: [T, NeuronTime, Xest] = fit_nnmf(X, 3, r);
%

% Canonical layout of data tensor X is [neurons x time x trials]
[num_neurons, num_samples, num_trials] = size(X);

switch unfold_dim
    case 1 % Neuron unfolding. (Fast-axis=Samples; Slow-axis=Trials)
        fprintf('Applying NEURAL unfolding; Xu=[num_neurons x num_samples(fast)*num_trials(slow)]\n');
        unfold_p = [1 2 3];
    case 3 % Trial unfolding. (Fast-axis=Samples, Slow-axis=Neurons)
        fprintf('Applying TRIAL unfolding; Xu=[num_trials x num_samples(fast)*num_neurons(slow)]\n');
        unfold_p = [3 2 1];
    otherwise
        error('fit_nnmf: Cannot unfold X along dim=%d!', unfold_dim);
end

% Generate the unfolded matrix along the requested direction
Xu = permute(X, unfold_p);
[~, num_fast, num_slow] = size(Xu);
Xu = reshape(Xu, [size(Xu,1), num_fast*num_slow]);

% Perform NMF on the unfolded matrix
options = statset('nnmf');
options.Display = 'iter'; % Verbose output
options.MaxIter = 200;
[A,B] = nnmf(Xu, r, 'options', options,...
             'algorithm', 'als');
                
Xest = A*B; % Reconstructed tensor

% Invert the permutation of dimensions of X back into canonical form
Xest = ipermute(Xest, unfold_p);
Xest = reshape(Xest, [num_neurons, num_samples, num_trials]);

B = reshape(B, [r, num_fast, num_slow]);
B = permute(B, [2 3 1]); % [Fast-dim x Slow-dim x rank]
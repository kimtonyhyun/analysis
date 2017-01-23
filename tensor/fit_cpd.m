function [models, best_models, S] = fit_cpd(X,varargin)
% FIT_CPD, fits a series of cp decompositions of increasing rank and
% plots the percentage of variance explained as function of model
% complexity
%
%     [models,rsq] = fit_cpd(X,[num_starts,10],[min_rank,1],[max_rank,15])

% parse optional inputs
p = inputParser;
p.addParameter('num_starts', 10);
p.addParameter('min_rank', 15);
p.addParameter('max_rank', 15);
p.addParameter('verbose', true);
p.addParameter('method','cp_als');
p.addParameter('num_samples',@(r) ceil(10*r*log(r)));
p.parse(varargin{:});

ns = p.Results.num_starts;
min_rank = p.Results.min_rank;
max_rank = p.Results.max_rank;

% create tensor object
Xt = tensor(X);
normX = norm(Xt);

% create struct array
models(ns, max_rank) = struct('error', NaN, 'decomp', []);
for s = 1:ns
    for r = 1:max_rank
        models(s,r).error = NaN;
    end
end

% main loop
% iterate in random order for a better waitbar
for s = 1:ns
    if p.Results.verbose
        fprintf(['\nOPTIMIZATION RUN ' num2str(s) '\n\t fitting models, rank: '])
    end
    for r = min_rank:max_rank
        fprintf([num2str(r) '.'])
        switch p.Results.method
        case 'cp_nnals'
            decomp = normalize(cp_nnals(Xt, r, 'printitn', false));
        case 'cp_als'
            decomp = normalize(cp_als(Xt, r, 'printitn', false));
        case 'cprand'
            ns = p.Results.num_samples(r);
            decomp = normalize(cprand(Xt, r, 'printitn', false, 'num_samples', ns, 'fft', 1));
        otherwise
            error('fitting method not recognized');
        end

    	models(s,r).error = sqrt(normX^2 + norm(decomp)^2 - 2*innerprod(Xt, decomp)) / normX;
    	models(s,r).decomp = decomp;
    end
end

if p.Results.verbose
    fprintf('\n')
end

best_models(1, max_rank) = struct('error', NaN, 'decomp', []);
for r = 1:(min_rank-1)
    best_models(r).error = NaN;
end
for r = min_rank:max_rank
    [~,s] = min([models(:,r).error]);
    best_models(r) = models(s,r);
end


S = nan(ns, ns, max_rank);
disp('Calculating pairwise similarities between model fits...')
for r = min_rank:max_rank
    disp(['Rank ', num2str(r) ' models.'])
    S(:,:,r) = cpd_pairwise_similarities(models(:,r));
end


function S = cpd_pairwise_similarities(models)
% CPD_PAIRWISE_SIMILARITIES, computes a similarity score for all pairs of
% fitted cpd models.
%
%     [similarity_matrix] = cpd_pairwise_similarities(models)

n_models = length(models);
S = zeros(n_models);

n_iter = n_models + (n_models^2)/2;
iter = 0;
prog = 0.1;

for a = 1:length(models)
    for b = a:length(models)
        S(a,b) = score(models(a).decomp, models(b).decomp, 'greedy', greedy);
        S(b,a) = S(a,b);
        
        % display progress
        iter = iter+1;
        if iter/n_iter > prog
        	disp([num2str(prog*100) '% done.'])
        	prog = prog+0.1;
        end
    end
end

disp('100% done.')

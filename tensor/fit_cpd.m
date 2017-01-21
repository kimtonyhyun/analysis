function [models, best_models] = fit_cpd(X,varargin)
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

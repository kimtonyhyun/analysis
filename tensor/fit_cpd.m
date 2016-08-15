function [cpd,rsq] = fit_cpd(X,max_rank)
% CPD_SCREE, fits a series of cp decompositions of increasing rank and
% plots the percentage of variance explained as function of model
% complexity
%
%     [cpd_list,rsq] = fit_cpd(X,max_rank)

% TODO: make these specifiable arguments
ns = 10; % number of starts per rank
min_rank = 1;

% create tensor object
Xt = tensor(X);

% create struct array
nr = max_rank-min_rank+1;
cpd(ns*nr) = struct('rank',0,'Rsq',0,'factors',struct([]),'lambdas',[],'decomp',[]);

% main loop
idx = 1;
rsq = zeros(ns*nr,1);
rnk = zeros(ns*nr,1);
for r = min_rank:max_rank
	for s = 1:ns
        h = waitbar((idx-1)/(nr*ns));
		decomp = normalize(cp_als(Xt,r,'printitn',0));
		cpd(idx).rank = r;
		cpd(idx).Rsq = 1 - norm(Xt-full(decomp))/norm(Xt - mean(X(:)));
		cpd(idx).factors = struct('neuron',decomp.u{1},'time',decomp.u{2},'trial',decomp.u{3});
		cpd(idx).lambda = decomp.lambda;
        cpd(idx).decomp = decomp;
        rsq(idx) = cpd(idx).Rsq;
        rnk(idx) = r;
        idx = idx+1;
	end
end
close(h)

% return a ROW vector for easy iteration
cpd = cpd';

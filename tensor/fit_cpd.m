function [cpd] = fit_cpd(X)
% FIT_CPD, fits canonical polyadic decomposition (cpd) to a data tensor
% X, returning results in a struct.
%
%     cpd = FIT_CPD(X)

% TODO: make these specifiable arguments
ns = 10; % number of starts per rank
min_rank = 1;
max_rank = 10;
seed = 1234; % for reproducible results

% dimensions
[nn,nt,nk] = size(X);

% create tensor object
Xt = tensor(X);

% create struct array
nr = max_rank-min_rank+1;
cpd(ns,nr) = struct('rank',0,'Rsq',0,'factors',struct([]),'lambdas',[],'decomp',[])

% main loop
ri = 1;
for r = min_rank:max_rank
	ri = r-min_rank+1;
	for s = 1:ns
		decomp = normalize(cp_als(Xt,r));
		cpd(s,ri).rank = r;
		cpd(s,ri).Rsq = 1 - norm(Xt-full(decomp))/norm(Xt - mean(X(:)));
		cpd(s,ri).factors = struct('neuron',decomp.u{1},'time',decomp.u{2},'trial',decomp.u{3});
		cpd(s,ri).lambda = decomp.lambda;
        cpd(s,ri).decomp = decomp;
	end
end


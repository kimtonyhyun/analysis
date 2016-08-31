function [cpd,rsq] = fit_cpd(X,varargin)
% FIT_CPD, fits a series of cp decompositions of increasing rank and
% plots the percentage of variance explained as function of model
% complexity
%
%     [cpd_list,rsq] = fit_cpd(X,[num_starts,10],[min_rank,1],[max_rank,15])

% parse optional inputs
p = inputParser;
p.addParameter('num_starts', 10);
p.addParameter('min_rank', 1);
p.addParameter('max_rank', 15);
p.addParameter('printitn',0);
p.parse(varargin{:});

ns = p.Results.num_starts;
min_rank = p.Results.min_rank;
max_rank = p.Results.max_rank;

% create tensor object
Xt = tensor(X);

% create struct array
nr = max_rank-min_rank+1;
cpd(1,ns*nr) = struct('rank',0,'Rsq',0,'factors',struct([]),'lambdas',[],'decomp',[]);

% main loop
rsq = zeros(ns*nr,1);

% iterate in random order for a better waitbar
itercount = 0;
for a = randperm(ns*nr)
    h = waitbar(itercount/(nr*ns));
    [r,~] = ind2sub([nr,ns],a);
	decomp = normalize(cp_als(Xt,min_rank+r-1,'printitn',p.Results.printitn));

	cpd(a).rank = min_rank+r-1;
	cpd(a).Rsq = 1 - norm(Xt-full(decomp))/norm(Xt - mean(X(:)));
	cpd(a).factors = struct('neuron',decomp.u{1},'time',decomp.u{2},'trial',decomp.u{3});
	cpd(a).lambda = decomp.lambda;
    cpd(a).decomp = decomp;
    rsq(a) = cpd(a).Rsq;
    itercount = itercount+1;
end
close(h)

% return a ROW vector for easy iteration
cpd = cpd';

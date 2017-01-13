function cpd = fit_cpd(X,varargin)
% FIT_CPD, fits a series of cp decompositions of increasing rank and
% plots the percentage of variance explained as function of model
% complexity
%
%     [cpd_list,rsq] = fit_cpd(X,[num_starts,10],[min_rank,1],[max_rank,15])

% parse optional inputs
p = inputParser;
p.addParameter('num_starts', 10);
p.addParameter('min_rank', 15);
p.addParameter('max_rank', 15);
p.addParameter('printitn',true);
p.addParameter('method','cp_als');
p.parse(varargin{:});

ns = p.Results.num_starts;
min_rank = p.Results.min_rank;
max_rank = p.Results.max_rank;

% create tensor object
Xt = tensor(X);
normX = norm(Xt);

% create struct array
nr = max_rank-min_rank+1;
cpd(max_rank, ns) = struct('error',0,'decomp',[]);

% main loop
% iterate in random order for a better waitbar
itercount = 0;
for s = 1:ns
    for r = min_rank:max_rank
        h = waitbar(itercount/(nr*ns));
        
        switch p.Results.method
        case 'cp_nnals'
            decomp = normalize(cp_nnals(Xt,min_rank+r-1,'printitn',p.Results.printitn));
        case 'cp_als'
            decomp = normalize(cp_als(Xt,min_rank+r-1,'printitn',p.Results.printitn));
        case 'cprand'
            decomp = normalize(cprand(Xt,min_rank+r-1,'printitn',p.Results.printitn));
        end

    	cpd(s,r).error = sqrt(normX^2 + norm(decomp)^2 - 2*innerprod(Xt, decomp));
    	cpd(s,r).decomp = decomp;
    end
end
close(h)

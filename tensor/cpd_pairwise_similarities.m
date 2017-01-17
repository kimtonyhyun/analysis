function S = cpd_pairwise_similarities(models, varargin)
% CPD_PAIRWISE_SIMILARITIES, computes a similarity score for all pairs of
% fitted cpd models.
%
%     [similarity_matrix] = cpd_pairwise_similarities(models)

n_models = length(models);
S = zeros(n_models);

disp('Calculating pairwise similarities...')

n_iter = n_models + (n_models^2)/2;
iter = 0;
prog = 0.1;

for a = 1:length(models)
    for b = a:length(models)
        S(a,b) = score(models(a).decomp, models(b).decomp, varargin{:});
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

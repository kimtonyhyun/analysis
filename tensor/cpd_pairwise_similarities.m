function S = cpd_pairwise_similarities(models, varargin)
% CPD_PAIRWISE_SIMILARITIES, computes a similarity score for all pairs of
% fitted cpd models.
%
%     [similarity_matrix] = cpd_pairwise_similarities(models)

S = zeros(length(models));

fprintf('Calculating pairwise similarities...\n')
for a = 1:length(models)
    fprintf(repmat(' ', 1, a))
    for b = a:length(models)
        fprintf('*')
        S(a,b) = score(models(a).decomp, models(b).decomp, varargin{:});
        S(b,a) = S(a,b);
    end
    fprintf('\n')
end


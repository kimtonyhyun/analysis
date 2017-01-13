function visualize_rank(cpd_list)
% VISUALIZE_RANK, plots a scree plot given a vector of cpd fits
%
%     visualize_rank(cpd_list)

% collect rank and rsq of each model
N = length(cpd_list);
ranks = zeros(N,1);
rsq = zeros(N,1);
for i = 1:N
    ranks(i) = size(cpd_list(i).decomp.u{1}, 2);
    rsq(i) = cpd_list(i).error;
end
unique_rank = unique(ranks);
n_unique = length(unique_rank);

% average rsq for each rank
avg_rsq = zeros(n_unique,1);
for i = 1:n_unique
    avg_rsq(i) = mean(rsq(ranks == unique_rank(i)));
end

% make figure
figure(); hold on
plot(unique_rank, 1-avg_rsq, '-r', 'linewidth', 2)
plot(ranks, 1-rsq, '.k', 'markersize', 20)
xlabel('rank of model')
ylabel('unexplained variance')
ylim([0,1])


% iterate over ranks
figure(); hold on
for i = 1:n_unique
    % identify models with this rank
    r = unique_rank(i);
    idx = find(ranks == r);
    
    % enumerate all pairs of models
    cidx = combnk(idx, 2);
    
    % calculate pairwise scores
    n_comb = size(cidx, 1);
    sc = zeros(n_comb, 1);
    for j = 1:n_comb
        model_1 = cpd_list(cidx(j,1)).decomp;
        model_2 = cpd_list(cidx(j,2)).decomp;
        sc(j) = score(model_1, model_2, 'greedy', true);
    end
    
    % plot
    jit = randn(n_comb, 1)*0.1;
    plot(ones(n_comb,1)*r + jit, sc, 'ok')
end

xlabel('model rank')
ylabel('model similarity')

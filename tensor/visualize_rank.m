function visualize_rank(cpd_list)
% VISUALIZE_RANK, plots a scree plot given a vector of cpd fits
%
%     visualize_rank(cpd_list)

% collect rank and rsq of each model
N = length(cpd_list);
rnk = zeros(N,1);
rsq = zeros(N,1);
for i = 1:N
    rnk(i) = cpd_list(i).rank;
    rsq(i) = cpd_list(i).Rsq;
end

% average rsq for each rank
rnk_uniq = unique(rnk)';
n = length(rnk_uniq);
avg_rsq = zeros(n,1);
for i = 1:n
    avg_rsq(i) = mean(rsq(rnk == rnk_uniq(i)));
end

% make figure
figure(); hold on
plot(unique(rnk),1-avg_rsq,'-r','linewidth',2)
plot(rnk,1-rsq,'.k','markersize',20)
xlabel('rank of model')
ylabel('unexplained variance')
ylim([0,1])

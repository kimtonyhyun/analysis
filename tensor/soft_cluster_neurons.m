function [cluster_ids, perm] = soft_cluster_neurons(factor)
% Returns soft-clustering of data based on CP decomposition results
%
% Input:
%   factor: N x R matrix of nonnegative data
%       Each neuron is a row; features in columns
%
% Output:
%   cluster_ids
%   perm: Permutation / ordering of rows by soft clustering

% Cluster based on score of maximum absolute value
max_num_clusters = size(factor, 2);
[scores, cluster_ids] = max(abs(factor), [], 2);
 
% Re-sort within each cluster
perm = [];
for k = 1:max_num_clusters
    inds_in_cluster = find(cluster_ids == k);
    scores_in_cluster = scores(inds_in_cluster);
    [~, ind] = sort(scores_in_cluster, 'descend');
    perm = [perm; inds_in_cluster(ind)]; %#ok<AGROW>
end
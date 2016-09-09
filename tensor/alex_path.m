[X, neuron_map, trial_map] = export_multiday_traces(md,'start','west','end','south');
daylist = trial_map(:,1);

X = timewarp(X);
for c = 1:size(X,1)
    for day = unique(daylist)'
        k = daylist == day;
        x = X(c,:,k);
        X(c,:,k) = x ./ (max([1, max(abs(x(:)))]));
    end
end

% find some rank 15 models
[cpd_list,rsq] = fit_cpd(X);
save('alex_path_cpdlist','cpd_list');

% % make a scree plot
% visualize_rank(cpd_list);

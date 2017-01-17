function [ln, sc] = cpd_scree_plot(models)
% CPD_SCREE_PLOT, plots a scree plot given struct array of model fits
%
%     models = fit_cpd(data, ...)
%     CPD_SCREE_PLOT(MODELS)
%

% collect rank and rsq of each model
n_replicates = size(models, 1);
max_rank = size(models, 2);

% scree plot
figure();
hold on
sc = [];
for r = 1:max_rank
    err = [models(:,r).error];
    mean_err(r) = mean(err); %#ok<AGROW>
    sc_ = plot(r*ones(n_replicates,1), err, '.k', 'markersize', 20);
    sc = [sc; sc_]; %#ok<AGROW>
end
rnks = find(~isnan(mean_err));
ln = plot(rnks, mean_err(rnks), '-', 'linewidth', 2, 'color', [1.0 0.4 0.4]);
uistack(ln, 'bottom')
xlabel('rank of model')
ylabel('normalized error')
ylim([0,1])

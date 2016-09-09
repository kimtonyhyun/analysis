function [B,cv_rate] = visualize_cpd_regress(cpd, meta)
% CPD_PREDICT(cpd, meta, trial_type)
%
% Uses logistic regression to predict trial metadata from cpd trial factors.
%
% model = cpd_predict(cpd, meta, 'start')
% model = cpd_predict(cpd, meta, 'end')
% model = cpd_predict(cpd, meta, 'correct')

% number of factors
nf = length(cpd.lambda);

% dependent variables to visualize
dep_vars = {'start', 'end', 'correct'};
ndv = length(dep_vars);

% storage for regresion coeffs and cv_rate
B = zeros(nf+1,ndv);
cv_rate = zeros(1,ndv);

for a = 1:ndv
    [B(:,a), cv_rate(a)] = cpd_regress_trial(cpd, meta, dep_vars{a});
end

% make plot showing strength of regression coefficients
figure()
bar(0:15,abs(B))
set(gca,'Xtick',0:15)
legend({ sprintf('start location (%1.3f accuracy)',cv_rate(1)), ...
         sprintf('end location (%1.3f accuracy)',cv_rate(2)),...
         sprintf('correct choice (%1.3f accuracy)',cv_rate(3)) });
xlim([-1 16])
ylabel('absolute value of regression coefficients')
xlabel('trial factors (0 = intercept term)')

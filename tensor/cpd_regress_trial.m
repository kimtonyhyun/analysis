function [B,cv_rate] = cpd_regress_trial(cpd, meta, trial_type)
% CPD_PREDICT(cpd, meta, trial_type)
%
% Uses logistic regression to predict trial metadata from cpd trial factors.
%
% model = cpd_predict(cpd, meta, 'start')
% model = cpd_predict(cpd, meta, 'end')
% model = cpd_predict(cpd, meta, 'correct')

% independent, dependent vars
X = cpd.factors.trial;
metadata = meta.(trial_type);
labels = unique(metadata);
Y = zeros(length(metadata),1);

for i = 1:length(metadata)
    for j = 1:length(labels)
        if metadata{i} == labels{j}
            Y(i) = j-1;
        end
    end
end

% find leave-one-out error rate
cv_rate = mean(crossval(@fit_predict,X,Y,'leaveout',1));

% fit the full model, return coefficients
B = glmfit(X,Y,'binomial','link','logit');


function num_correct = fit_predict(xT,yT,xt,yt)
    % train on (xT,yT) then predict yt using xt
    B = glmfit(xT,yT,'binomial','link','logit');
    yhat = glmval(B,xt,'logit');
    num_correct = sum(abs(yt-yhat)<0.01);

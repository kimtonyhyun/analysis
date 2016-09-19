function [B,cv_rate] = cpd_regress_trial(cpd, meta, var)
% CPD_PREDICT(cpd, meta, trial_type)
%
% Uses multinomial logistic regression to predict day, strategy, start arm,
% etc. from trial factors.

if strcmp(var, 'strategy')
    % remove non-labelled trials
    idx = strcmp('allo-north',meta.strategy) | ...
          strcmp('allo-south',meta.strategy) | ...
          strcmp('ego-right',meta.strategy) | ...
          strcmp('ego-left',meta.strategy);

    X = cpd.factors.trial(idx,:);
    Y = categorical(meta.strategy(idx));
else
    X = cpd.factors.trial;
    Y = categorical(meta.(var));
end
    
% find leave-one-out error rate
cv_rate = mean(crossval(@fit_predict,X,Y,'kfold',10));

% fit the full model, return coefficients
B = mnrfit(X,Y);

function frac_correct = fit_predict(xT,yT,xt,yt)
    % train on (xT,yT) then predict yt using xt
    B = mnrfit(xT,yT);
    yhat = mnrval(B,xt);
    [~,yi] = max(yhat,[],2);
    yl = categories(yT);
    frac_correct = mean(yl(yi) == yt);

function [B,cv_rate] = cpd_regress_strat(cpd, meta)
% CPD_PREDICT(cpd, meta, trial_type)
%
% Uses multinomial logistic regression to predict strategy from trial factors.

% independent, dependent vars
idx = strcmp('allo-north',meta.strategy) | ...
      strcmp('allo-south',meta.strategy) | ...
      strcmp('ego-right',meta.strategy) | ...
      strcmp('ego-left',meta.strategy);

X = cpd.factors.trial(idx,:);
Y = categorical(meta.strategy(idx));

% find leave-one-out error rate
cv_rate = mean(crossval(@fit_predict,X,Y,'leaveout',1));

% fit the full model, return coefficients
B = mnrfit(X,Y);

function num_correct = fit_predict(xT,yT,xt,yt)
    % train on (xT,yT) then predict yt using xt
    B = mnrfit(xT,yT);
    yhat = mnrval(B,xt);
    [~,yi] = max(yhat,[],2);
    yl = categories(yT);
    num_correct = sum(yl(yi) == yt);

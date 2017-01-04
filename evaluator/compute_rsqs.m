function [Rsqs, Rsqs_n] = compute_rsqs(X, Xest)

assert(all(size(X)==size(Xest)),...
    'Raw and fitted data do not have the same dimensions!');

[num_neurons, ~, num_trials] = size(X);

Rsqs = zeros(num_neurons, num_trials);
Rsqs_n = zeros(num_neurons, 1);

for nid = 1:num_neurons
    for tid = 1:num_trials
        tr = X(nid,:,tid);
        tr_est = Xest(nid,:,tid);
        
        Rsqs(nid,tid) = compute_rsq(tr, tr_est);
    end
    
    tr = X(nid,:,:); tr = tr(:);
    tr_est = Xest(nid,:,:); tr_est = tr_est(:);
    
    Rsqs_n(nid) = compute_rsq(tr, tr_est);
end

end % compute_rsqs

function rsq = compute_rsq(y, y_est)

    SS_res = sum((y-y_est).^2);
    SS_tot = sum((y-mean(y)).^2);

    rsq = 1 - SS_res / SS_tot;

end % compute_rsq
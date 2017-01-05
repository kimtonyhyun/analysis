function scores = compute_all_scores(X, Xest)

% Definition of score struct
%--------------------------------------------------
score = struct('name', '',...       % String name of score
               'sort_order', '',... % To order from best to worst
               'vals', [],...      % Scores: [num_neurons x num_trials]
               'vals_n', []);      % Score per neuron: [num_neurons x 1]

num_score_types = 3;
scores = repmat(score, num_score_types, 1);
            
scores(1).name = 'Rsq';
[scores(1).vals, scores(1).vals_n] = compute_rsqs(X, Xest);
scores(1).sort_order = 'descend';

scores(2).name = 'Residual';
[scores(2).vals, scores(2).vals_n] = compute_residuals(X, Xest);
scores(2).sort_order = 'ascend';

scores(3).name = 'Peak-to-peak';
Xmin = squeeze(min(X,[],2));
Xmax = squeeze(max(X,[],2));
scores(3).vals = Xmax - Xmin;
scores(3).vals_n = max(Xmax,[],2) - min(Xmin,[],2);
scores(3).sort_order = 'descend';
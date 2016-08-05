clear; close all; clc;

if ismac
    p = '/Users/alex/Dropbox/strategy_switch/data/c11m1d12';
elseif isunix
    p = '/home/alex/Dropbox/strategy_switch/data/c11m1d12';
else
    error('Platform not supported')
end

% do tensor analysis
sources = gen_data_sources(p);
session = DaySummary(sources, [p,'/cm01']);
[X,cell_idx,trial_idx] = session_tensor(session);
X = soft_normalize(X);
cpd = fit_cpd(X);
% factor_plots(session,cpd(10,10),trial_idx);
scree_plot(cpd);


% Xest = full(cpd(10,10).decomp);
% neuron_fit_plots(session,X,Xest.data,cell_idx,trial_idx);
neuron_factor_plots(session,cpd(10,10),trial_idx,'start')

% return to original directory
cd(p)

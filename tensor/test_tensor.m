clear; close all; clc;

% change to data folder
p = pwd();
cd('/Users/alex/Dropbox/strategy_switch/data/c11m1d12')

% do tensor analysis
sources = gen_data_sources;
session = DaySummary(sources, 'cm01');
[X,cell_idx,trial_idx] = session_tensor(session);
X = soft_normalize(X);
X = X - mean(X(:));
cpd = fit_cpd(X);
factor_plots(session,cpd(10,10),trial_idx);
scree_plot(cpd);

Xest = full(cpd(10,10).decomp);
neuron_fit_plots(session,X,Xest.data,cell_idx,trial_idx);

% return to original directory
cd(p)

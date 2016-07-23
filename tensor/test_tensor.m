clear; close all; clc;

cd('/home/alex/data/strategy_switch/c11m1d12')

sources = gen_data_sources;
session = DaySummary(sources, 'cm01');
[X,cell_idx,trial_idx] = session_tensor(session);
X = soft_normalize(X);
cpd = fit_cpd(X);
factor_plots(session,cpd(10,10),trial_idx);
scree_plot(cpd);

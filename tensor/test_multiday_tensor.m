clear; close all; clc;

if ismac
    p_m1d12 = '/Users/alex/Dropbox/strategy_switch/data/c11m1d12';
    p_m1d13 = '/Users/alex/Dropbox/strategy_switch/data/c11m1d13';
    p_m1d14 = '/Users/alex/Dropbox/strategy_switch/data/c11m1d14';
elseif isunix
    p_m1d12 = '/home/alex/Dropbox/strategy_switch/data/c11m1d12';
    p_m1d13 = '/home/alex/Dropbox/strategy_switch/data/c11m1d13';
    p_m1d14 = '/home/alex/Dropbox/strategy_switch/data/c11m1d14';
else
    error('Platform not supported')
end

src_m1d12 = gen_data_sources(p_m1d12);
src_m1d13 = gen_data_sources(p_m1d13);
src_m1d14 = gen_data_sources(p_m1d14);

m1d12 = DaySummary(src_m1d12, [p_m1d12,'/cm01']);
m1d13 = DaySummary(src_m1d13, [p_m1d13,'/cm01']);
m1d14 = DaySummary(src_m1d14, [p_m1d14,'/cm01']);
ds_list = {12, m1d12; 13, m1d13; 14, m1d14};

[X,cell_idx,trial_idx] = multiday_tensor(ds_list);

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

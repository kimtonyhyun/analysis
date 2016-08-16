%% This is a hacky script to test things ... Tony disapproves

clear; close all; clc;

% replace with your path, or comment this out and load data yourself
if ismac
    p_m1d12 = '/Users/alex/Dropbox/strategy_switch/data/c11m1d12';
    p_m1d13 = '/Users/alex/Dropbox/strategy_switch/data/c11m1d13';
    p_m1d14 = '/Users/alex/Dropbox/strategy_switch/data/c11m1d14';
    figdir = '/Users/alex/Dropbox/strategy_switch/figs/';
    load /Users/alex/Dropbox/strategy_switch/data/match_m1d12-14.mat
elseif isunix
    p_m1d12 = '/home/alex/Dropbox/strategy_switch/data/c11m1d12';
    p_m1d13 = '/home/alex/Dropbox/strategy_switch/data/c11m1d13';
    p_m1d14 = '/home/alex/Dropbox/strategy_switch/data/c11m1d14';
    figdir = '/home/alex/Dropbox/strategy_switch/figs/';
    load /home/alex/Dropbox/strategy_switch/data/match_m1d12-14.mat
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

match_list = {12, 13, m_12to13, m_13to12;
              12, 14, m_12to14, m_14to12;
              13, 14, m_13to14, m_14to13};

md = MultiDay(ds_list,match_list);

%% The real tensor stuff happens here
md = MultiDay({12,m1d12},{});
[X, neuron_map, trial_map] = export_multiday_traces(md,'all');

% TODO: standardize before warp?
% X = standardize(X,trial_map);

% preprocessing.
X = timewarp(X);
for c = 1:size(X,1)
    x = X(c,:,:);
    X(c,:,:) = x ./ (max([1, max(abs(x(:)))]));
end

% make a scree plot (up to rank 10)
[cpd_list,rsq] = fit_cpd(X,10);
scree_cpd(cpd_list);

% pick best cpd to analyze further
[~,i] = max(rsq);
cpd = cpd_list(i);
Xest = full(cpd.decomp);
Xest = Xest.data;

% plot single-figure summary of all factors
neuron_factor_plots(cpd,md,trial_map)

% plot fit across neurons
outdir = [figdir, 'm1/neuron_viz/'];
visualize_fit(X,Xest,1,md,trial_map,outdir);

% plot fit across trials
outdir = [figdir, 'm1/trial_viz/'];
visualize_fit(X,Xest,3,md,trial_map,outdir);

% plot residuals
visualize_resids(X,Xest,md,trial_map);


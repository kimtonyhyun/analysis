function visualize_neuron_ktensor(decomp, meta, trialcolor)
% VISUALIZE_NEURON_KTENSOR, wraps VISUALIZE_KTENSOR for multi-day data
%
%     [Ax, BigAx, FigHandle] = VISUALIZE_NEURON_KTENSOR(X)
%
%     Parameters
%     ----------
%     trialcolor : 'start', 'error'

% color for the trials
c = cell(1,3);
if nargin >=3
    c{3} = get_trial_colors(meta.(trialcolor));
end

% the neuron order is arbitrary, so permute/sort it along the first factor
prm = [true false false];

% bar plot for neuron factors
% line plot for within trial factors
% line and scatter plot for across trial factors
plt = {'bar', 'line', 'scatter'};

% plot titles
nm = {'neuron factors', 'time factors', 'trial factors'};

% formatting
lspc = {[], '-r', '-k'};
lw = [0 2 1];

visualize_ktensor(decomp, 'c', c, 'plots', plt, ...
                  'permute', prm, 'linewidth', lw, ...
                  'names', nm, 'linespec', lspc);

function colors = get_trial_colors(labels)
% trial coloring labels

% convert labels to categorical array
if ~iscategorical(labels)
    labels = categorical(labels);
end
cats = categories(labels);
nc = length(cats);

% generate colormap
gr = 0.68;%618033988749895;
h = 0;
cm = zeros(nc,3);
for c = 1:length(cats)
    cm(c,:) = hsv2rgb(h,1,1);
    h = mod(h+gr,1);
end

% assign colors to datapoints
K = length(labels);
colors = zeros(K,3);
for k = 1:K
    colors(k,:) = cm(cats == labels(k),:);
end

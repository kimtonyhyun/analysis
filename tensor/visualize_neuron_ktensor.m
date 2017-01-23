function [Ax, G, BigAx] = visualize_neuron_ktensor(decomp, meta, trialcolor)
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
    c{3} = categorical_colors(meta.(trialcolor));
end

% the neuron order is arbitrary, so permute/sort it along the first factor
[~, n_idx] = sort(sum(abs(decomp.u{1}), 2), 'descend');
decomp.u{1} = decomp.u{1}(n_idx, :) .* repmat(decomp.lambda', length(n_idx), 1);

% bar plot for neuron factors
% line plot for within trial factors
% line and scatter plot for across trial factors
plt = {'bar', 'line', 'scatter'};

% plot titles
nm = {'neuron factors', 'time factors', 'trial factors'};

% formatting
lspc = {[], '-r', '-k'};
lw = [0 2 1];

[Ax, G, BigAx] = visualize_ktensor(decomp, 'c', c, 'plots', plt, ...
                  'linewidth', lw, 'link_yax', [true, true, true], ...
                  'title', nm, 'linespec', lspc);


set(Ax(1:end-1,:), 'YTickLabel', [])
set(Ax, 'YAxisLocation', 'left')
set(Ax(end, :), 'TickDir', 'out')

function neuron_factor_plots(decomp,trial_idx,varargin)

% tensor dimensions (neurons x time x trial)
factors = decomp.factors;
nn = size(factors.neuron,1);
nr = size(factors.neuron,2);
nk = length(trial_idx);

% plot factors in order of decreasing variability across trials
peak_idx = zeros(nr,1);
for r = 1:nr
    peak_idx(r) = std(factors.trial(:,r));
end
[~,fo] = sort(peak_idx,'descend');

% color the trials
if nargin == 3
    colormode = 'none';
else
    colormode = varargin{1};
end

trial_colors = zeros(nk,3);
if ~strcmp(colormode,'none')
    [i1,i2,labels] = ten_filter_trials(session,trial_idx,colormode);
    for k = transpose(i1)
        trial_colors(k,:) = [0 0 1];
    end
    for k = transpose(i2)
        trial_colors(k,:) = [1 0 0];
    end
end

% plot trials by true order
if nargin < 5
    axmode = 'order';
else
    axmode = varargin{2};
end

switch axmode
    case 'order'
        trial_ax = 1:nk;
    case 'number'
        trial_ax = trial_idx;
    otherwise
        error('unsupported trial axis')
end

% make the figure
figure()
subplot(1,3,1)
plotmatrix(factors.neuron(:,fo))
title('neuron factors')

subplot(1,3,2)
set(gca,'Visible','off')
pos = get(gca,'Position');
width = pos(3);
height = pos(4)/10;
space = .02; % 2 percent space between axes
pos(1:2) = pos(1:2) + space*[width height];

ax = gobjects(10);
yl = 1.01*max(abs(factors.trial(:)));
for r = 1:nr
    axPos = [pos(1) pos(2)+(10-r)*height ...
                width*(1-space) height*(1-space)];
    ax(r) = axes('Position',axPos);
    hold on
    plot(trial_ax,factors.trial(:,fo(r)),'-k')
    scatter(trial_ax,factors.trial(:,fo(r)),20,trial_colors,'filled')
    set(gca,'xtick',[])
    ylim([-yl yl])
end

subplot(1,3,3)
set(gca,'Visible','off')
pos = get(gca,'Position');
width = pos(3);
height = pos(4)/10;
space = .02; % 2 percent space between axes
pos(1:2) = pos(1:2) + space*[width height];

ax = gobjects(10);
for r = 1:nr
    axPos = [pos(1) pos(2)+(10-r)*height ...
                width*(1-space) height*(1-space)];
    ax(r) = axes('Position',axPos);
    hold on
    plot(factors.time(:,fo(r)),'-k','linewidth',2)
    set(gca,'xtick',[])
end
title(ax(1),'time factors')
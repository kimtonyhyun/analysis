function cpd_factor_plots(cpd,md,trial_map,varargin)

% parse optional inputs
params = inputParser;
params.addParameter('trialcolor', 'start', ...
                    @(x) any(validatestring(x,['start','error','none'])));
params.addParameter('trialax', 'order', ...
                    @(x) any(validatestring(x,['order','number'])));
params.addParameter('factor_order', 'lambda', ...
                    @(x) any(validatestring(x,['lambda','trialvar'])));
params.addParameter('neuron_plot', 'bars', ...
                    @(x) any(validatestring(x,['bars','plotmatrix'])));
params.addParameter('n_factors', cpd.rank);
params.addParameter('space', 0.1);
params.parse(varargin{:});
res = params.Results;
nf = res.n_factors;
space = res.space;

% tensor dimensions (neurons x time x trial)
factors = cpd.factors;
nn = size(factors.neuron,1);
nr = size(factors.neuron,2);
nk = size(trial_map,1);

% plot factors in order of decreasing variability across trials
switch res.factor_order
    case 'lambda'
        [~,fo] = sort(cpd.lambda,'descend');
    case 'trialvar'
        factvar = zeros(nr,1);
        for r = 1:nr
            factvar(r) = std(factors.trial(:,r));
        end
        [~,fo] = sort(factvar,'descend');
end

trial_colors = get_trial_colors(md, trial_map, res.trialcolor);

% plot trials by order or by true number
switch res.trialax
    case 'order'
        trial_ax = 1:nk;
    case 'number'
        trial_ax = trial_map(:,2);
        for d = sort(md.valid_days,'ascend')
            n = length(md.day(d).trials);
            idx = trial_map(:,1) > d;
            trial_ax(idx) = trial_ax(idx) + n;
        end
end

% make the figure
figure()

subplot(1,3,1);
switch res.neuron_plot
    case 'bars'
        [~,no] = sort(cpd.factors.neuron(:,fo(1)),'descend');
        set(gca,'Visible','off')
        pos = get(gca,'Position');
        width = pos(3);
        height = pos(4)/nf;
        pos(1:2) = pos(1:2) + space*[width height];

        ax = gobjects(nf);
        yl = 1.01*max(abs(factors.trial(:)));
        for r = 1:nf
            axPos = [pos(1) pos(2)+(nf-r)*height ...
                        width*(1-space) height*(1-space)];
            ax(r) = axes('Position',axPos);
            hold on
            bar(1:nn,factors.neuron(no,r))
            set(gca,'xtick',[],'xlim',([0,nn+1]),...
                    'ytick',[-0.1,0.1],'ylim',[-0.15,0.15])
        end
        title(ax(1),'neuron factors')
        
    case 'plotmatrix'
        plotmatrix(factors.neuron(:,fo))
        title('neuron factors')
end

subplot(1,3,2);
set(gca,'Visible','off')
pos = get(gca,'Position');
width = pos(3);
height = pos(4)/nf;
pos(1:2) = pos(1:2) + space*[width height];

ax = gobjects(nf);
yl = 1.01*max(abs(factors.trial(:)));
for r = 1:nf
    axPos = [pos(1) pos(2)+(nf-r)*height ...
                width*(1-space) height*(1-space)];
    ax(r) = axes('Position',axPos);
    hold on
    plot(trial_ax,factors.trial(:,fo(r)),'-k')
    scatter(trial_ax,factors.trial(:,fo(r)),20,trial_colors,'filled')
    set(gca,'xtick',[])
    ylim([-yl yl])
end
title(ax(1),'trial factors')

subplot(1,3,3); title('time factors')
set(gca,'Visible','off')
pos = get(gca,'Position');
width = pos(3);
height = pos(4)/nf;
pos(1:2) = pos(1:2) + space*[width height];

ax = gobjects(nf);
for r = 1:nf
    axPos = [pos(1) pos(2)+(nf-r)*height ...
                width*(1-space) height*(1-space)];
    ax(r) = axes('Position',axPos);
    plot(factors.time(:,fo(r)),'-k','linewidth',2)
    axis tight
    set(gca,'xtick',[])
    if mod(r,2) == 0
        set(gca,'YAxisLocation','right')
    end
    yl = get(gca,'ylim');
    ryl = round(yl,2);
    if ryl(2) > yl(2)
        yl(2) = ryl(2);
    end
    if ryl(1) < yl(1)
        yl(1) = ryl(1);
    end
    set(gca,'ytick',[ryl(1),ryl(2)],'ylim',yl)
end
title(ax(1),'time factors')
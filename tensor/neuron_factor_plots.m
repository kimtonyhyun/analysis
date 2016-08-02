function neuron_factor_plots(session,decomp,trial_idx)

% tensor dimensions (neurons x time x trial)
factors = decomp.factors;
nn = size(factors.neuron,1);
nr = size(factors.neuron,2);

[ie,iw] = ie_iw(session,trial_idx);

% peak_idx = zeros(nr,1)
% for r = 1:nr
%     [~,ir] = max(abs(factors.time(:,r)))
%     peak_idx(r) = ir
% end
% [~,fact_order] = sort(peak_idx)
% plotmatrix(factors.neuron(:,fact_order))

peak_idx = zeros(nr,1);
for r = 1:nr
    peak_idx(r) = std(factors.trial(:,r));
end
[~,fact_order] = sort(peak_idx);


figure()
subplot(1,3,1)
plotmatrix(factors.neuron(:,fact_order))
title('neuron factors')

subplot(1,3,2)
set(gca,'Visible','off')
pos = get(gca,'Position');
width = pos(3);
height = pos(4)/10;
space = .02; % 2 percent space between axes
pos(1:2) = pos(1:2) + space*[width height];

ax = gobjects(10);
yl = 0;
for r = 1:nr
    axPos = [pos(1) pos(2)+(10-r)*height ...
                width*(1-space) height*(1-space)];
    ax(r) = axes('Position',axPos);
    hold on
    plot(factors.trial(:,fact_order(r)),'-k')
    plot(ie,factors.trial(ie,fact_order(r)),'.r','markersize',20)
    plot(iw,factors.trial(iw,fact_order(r)),'.b','markersize',20)
    set(gca,'xtick',[])
    if yl < max(abs(ylim()))
        yl = max(abs(ylim()));
    end
end

for r = 1:nr
    set(ax(r),'ylim',[-yl yl])
end
title(ax(1),'trial factors')

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
    plot(factors.time(:,fact_order(r)),'-k','linewidth',2)
    set(gca,'xtick',[])
end
title(ax(1),'time factors')
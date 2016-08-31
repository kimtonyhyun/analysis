function visualize_factors(tnsrlist,md,trial_map,varargin)
% VISUALIZE_FACTORS, given a cell array of ktensors, plot the factors
%
%     [H, Ax, BigAx] = VISUALIZE_FACTORS(cpd)
%     [H, Ax, BigAx] = VISUALIZE_FACTORS(cpd, ['align', true])
%     [H, Ax, BigAx] = VISUALIZE_FACTORS(cpd, ['trialcolor', ['start','error','none']])
%     [H, Ax, BigAx] = VISUALIZE_FACTORS(cpd, ['order', ['order','number']])
%     [H, Ax, BigAx] = VISUALIZE_FACTORS(cpd, ['nfactors', cpd.rank])
%     [H, Ax, BigAx] = VISUALIZE_FACTORS(cpd, ['space', 0.1])
%     [H, Ax, BigAx] = VISUALIZE_FACTORS(cpd, ['figure', -1])

% check inputs
if ~iscell(tnsrlist) && ~isa(tnsrlist,'ktensor')
    tnsrlist
elseif ~iscell(tnsrlist)
    tnsrlist = {tnsrlist}
end

if ~isstruct(cpdlist)
    error('cpd must be a struct or struct array')
end

% parse optional inputs
params = inputParser;
params.addParameter('trialcolor', 'start', ...
                    @(x) any(validatestring(x,['start','error','none'])));
params.addParameter('trialax', 'order', ...
                    @(x) any(validatestring(x,['order','number'])));
params.addParameter('neuron_plot', 'bars', ...
                    @(x) any(validatestring(x,['bars','plotmatrix'])));
params.addParameter('nfactors', cpd.rank);
params.addParameter('align', true);
params.addParameter('space', 0.1);
params.addParameter('figure', -1);
params.parse(varargin{:});
res = params.Results;
nf = res.nfactors;

% get appropriate figure
if isnumeric(res.figure) && res.figure < 1
    fh = figure()
elseif isgraphics(res.figure)
    fh = figure(res.figure)
else
    error('invalid figure handle')
end

% tensor dimensions (neurons x time x trial)
factors = cpd.factors;
nn = size(factors.neuron,1);
nr = size(factors.neuron,2);
nk = size(trial_map,1);

% plot factors in order of decreasing variability across trials
[~,fo] = sort(cpd(1).lambda,'descend');

% determine color of trial factor plots
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
fh = figure()

[Ax,BigAx] = setup_axes(nf, res.space)

% if there are multiple fits, try to align factors
% and plot them on the same axes
for idx = 1:length(cpdlist)

    % cpd to plot
    cpd = cpdlist(idx)

    % align all cpds to the first one
    if idx > 1 && res.align
        % use score to align the CPDs
        [~,nd] = score(cpd(1).decomp, cpd(idx).decomp)
        factors = nd.factors
    else
        % don't do alignment
        factors = cpd.factors
    end

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



% LOCAL FUNCTIONS %
function [Ax,BigAx] = setup_axes(nfactors, space)

    % allocate storage
    Ax = gobjects(nfactors,3);
    BigAx = gobjects(3);
    
    % setup axes
    for bi = 1:3

        % invisible subplot bounding box
        BigAx[bi] = subplot(1,3,bi);
        set(BigAx[bi],'Visible','off')
        pos = get(BigAx[bi],'Position');
        w = pos(3);
        h = pos(4)/nf;
        pos(1:2) = pos(1:2) + space*[w h];

        % subaxes
        for si = 1:nfactors
            axPos = [pos(1) pos(2)+(nf-si)*h w*(1-s) h*(1-s)];
            Ax(si,bi) = axes('Position',axPos);
        end
    end
end



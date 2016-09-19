function visualize_fit(X,Xest,dim,md,trial_map,outdir)
% VISUALIZE_FIT, compare the cpd fit (Xest) to the real data (X) along mode
% dim.
%
%     VISUALIZE_FIT(X,Xest,dim,ds_list,trial_map)

figure()

if nargin==5
    outdir = '';
end

if dim == 1
    ystr = 'trial';
elseif dim == 3
    Xest = permute(Xest,[3 2 1]);
    X = permute(X,[3 2 1]);
    ystr = 'neuron';
else
    error('cant vizualize that')
end

s = randsample(1:size(X,3),25);
xlimits = [1, size(X,2)];

nframes = size(X,1);
for n = 1:nframes
    clf()
    ylimits = [Inf,-Inf];
    for k = 1:25
        
        % select day summary
        if dim == 1
            trialX = k;
        else
            trialX = n;
        end
        day = trial_map(trialX,1);
        trialDAY = trial_map(trialX,2);
        
        ax(k) = subplot(5,5,k);
        cla; hold on
        if strcmp(md.day(day).trials(trialDAY).start, 'east')
            plot(X(n,:,s(k)),'-','linewidth',3,'color',[0 0.7 1])
        elseif strcmp(md.day(day).trials(trialDAY).start, 'west')
            plot(X(n,:,s(k)),'-r','linewidth',3)
        else
            warn('probe trial plotted')
        end
        plot(Xest(n,:,s(k)),'-k','linewidth',2)
        set(gca,'xtick',[])
        yl = ylim();
        ylimits(1) = min([yl(1), ylimits(1)]);
        ylimits(2) = max([yl(2), ylimits(2)]);
        xlim(xlimits)
        ylabel([ystr,' ', num2str(s(k))])
    end
    linkaxes(ax)
    ylim(ylimits)
    if isempty(outdir)
        pause
    else
        saveas(gcf, [outdir,sprintf('%03d',n)], 'pdf')
    end
end


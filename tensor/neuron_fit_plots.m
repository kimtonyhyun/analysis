function neuron_fit_plots(ds,X,Xest,cell_idx,trial_idx)

[ie,~] = ie_iw(ds,trial_idx);
trials = 1:4:100;
ylimits = [min(X(:)) max(X(:))];

figure('Position',[0 0 1200 1000])
for n = 1:size(X,1)
    cn = cell_idx(n);
    for k = 1:25
        subplot(5,5,k);
        cla; hold on
        if any(ie == trials(k))
            plot(X(n,:,trials(k)),'-r','linewidth',3)
        else
            plot(X(n,:,trials(k)),'-b','linewidth',3)
        end
        plot(Xest(n,:,trials(k)),'-k','linewidth',2)
        set(gca,'xtick',[],'ytick',[])
        ylim(ylimits)
        ylabel(['trial ', num2str(trials(k))])
    end 
    pause
%     in = input('Enter `q` to quit, any other button to continue:  ','s');
%     if strcmp(in,'q')
%         break
%     end
end

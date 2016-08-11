function visualize_fit(X,Xest_tnsr,dim)
% VISUALIZE_FIT, compare the cpd fit (Xest) to the real data (X) along mode
% dim.
%
%     VISUALIZE_FIT(X,Xest,dim)

if dim == 1
    Xest = Xest_tnsr.data;
    ystr = 'trial';
elseif dim == 3
    Xest = permute(Xest_tnsr,[3 2 1]);
    Xest = Xest.data;
    ystr = 'neuron';
else
    error('cant vizualize that')
end

s = randsample(1:size(X,3),25);

for n = 1:size(X,1)
    for k = 1:25
        subplot(5,5,k);
        cla; hold on
        plot(X(n,:,s(k)),'-b','linewidth',3)
%         if any(ie == s(k))
%             plot(X(n,:,s(k)),'-r','linewidth',3)
%         else
%             plot(X(n,:,s(k)),'-b','linewidth',3)
%         end
        plot(Xest(n,:,s(k)),'-k','linewidth',2)
        set(gca,'xtick',[],'ytick',[])
%         ylim(ylimits)
        ylabel([ystr,' ', num2str(s(k))])
    end 
    pause
%     in = input('Enter `q` to quit, any other button to continue:  ','s');
%     if strcmp(in,'q')
%         break
%     end
end

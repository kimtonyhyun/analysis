function scree_plot(cpd)
% SCREE_PLOT, shows unexplained variance as a function of number of factors. 
%
%     FACTOR_PLOTS(session,decomp,trial_idx)

vc = cpd(:);
n = length(vc);
x = zeros(n,1);
y = zeros(n,1);
for a = 1:n
    x(a) = cpd(a).rank;
    y(a) = 1-cpd(a).Rsq;
end

xlev = unique(x)';
nlev = length(xlev);
for a = 1:nlev
    yy(a) = mean(y(x == xlev(a)));
end

figure(); hold on
plot(xlev,yy,'-r')
plot(x,y,'.k')
xlabel('rank of model')
ylabel('unexplained variance')

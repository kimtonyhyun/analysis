function y_lims = tight_plot(x, y)
% Wrapper around "plot(x,y)" that automatically sizes xlim and ylim.
plot(x, y);
xlim(x([1 end]));
y_lims = [min(y) max(y)];
y_lims = y_lims + 0.1*diff(y_lims)*[-1 1];
ylim(y_lims);
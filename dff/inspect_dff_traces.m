function [tr2_dff, info] = inspect_dff_traces(ds_proj, ds_ls, idx, fps, varargin)
% DFF traces are most straightforward to compute when using simple
% projection to the motion-corrected imaging data (with minimal additional
% processing). However, least-squares projection provides better safeguards
% against trace contamination.
%
% This basic visualization function compares DFF traces obtained by simple
% projection vs. least-squares, as a sanity check that the DFF scaling of
% the two methods are comparable. In addition, the baseline (F0) estimate
% is shown for both projection and least squares traces.
%
% Returns: Results of 'compute_dff_trace' for the least squares trace
%

% Default colors
c1 = [0 0.447 0.741];
c2 = [0.85 0.325 0.098];

% Custom "subplot" command that leaves less unusued space between panels
sp = @(m,n,p) subtightplot(m, n, p, 0.05, 0.05, 0.05); % Gap, Margin-X, Margin-Y

tr1 = ds_proj.get_trace(idx);
[tr1_dff, info] = compute_dff_trace(tr1, varargin{:});
b1 = info.baseline;
nu1 = calculate_noise_level_nu(tr1_dff, fps);

tr2 = ds_ls.get_trace(idx);
[tr2_dff, info] = compute_dff_trace(tr2, varargin{:});
b2 = info.baseline;
nu2 = calculate_noise_level_nu(tr2_dff, fps);

num_frames = length(tr1);
t = 1/fps*(0:num_frames-1);

ax1 = sp(3,1,1);
plot(t, tr1, 'Color', c1);
hold on;
plot(t, b1, 'k-', 'LineWidth', 2);
hold off;
ylabel('Simple projection');
title(sprintf('Cell %d: FPS=%.1f Hz', idx, fps));
grid on;
set(ax1, 'TickLength', [0 0]);

ax2 = sp(3,1,2);
plot(t, tr2, 'Color', c2);
hold on;
plot(t, b2, 'k-', 'LineWidth', 2);
hold off;
ylabel('Least squares');
grid on;
set(ax2, 'TickLength', [0 0]);

ax3 = sp(3,1,3);
plot(t, tr1_dff, 'Color', c1);
hold on;
plot(t, tr2_dff, '--', 'Color', c2);
plot(t([1 end]), [0 0], 'k-', 'LineWidth', 2);
hold off;
ylabel('DFF');
grid on;
set(ax3, 'TickLength', [0 0]);
title(sprintf('Noise level \\nu = %.1f %%/Hz^{0.5} (LS)', nu2));

linkaxes([ax1 ax2 ax3], 'x');
xlim(t([1 end]));
xlabel('Time (s)');
zoom xon;

end

function compare_traces(ds_proj, ds_ls, idx, fps)
% DFF traces are most straightforward to compute when using simple
% projection to the motion-corrected imaging data (with minimal additional
% processing). However, least-squares projection provides better safeguards
% against trace contamination.
%
% This basic visualization function compares DFF traces obtained by simple
% projection vs. least-squares, as a sanity check that the DFF scaling of
% the two methods are comparable.

c1 = [0 0.447 0.741];
c2 = [0.85 0.325 0.098];

percentile_window = 3e3;

tr1 = ds_proj.get_trace(idx);
[tr1_dff, info] = compute_dff_trace(tr1, 'window', percentile_window);
b1 = info.baseline;
tr2 = ds_ls.get_trace(idx);
[tr2_dff, info] = compute_dff_trace(tr2, 'window', percentile_window);
b2 = info.baseline;

num_frames = length(tr1);
t = 1/fps*(0:num_frames-1);

ax1 = subplot(311);
plot(t, tr1, 'Color', c1);
hold on;
plot(t, b1, 'k-');
hold off;
ylabel('Projection');
title(sprintf('Cell %d: FPS=%.1f Hz, Baseline window=%d frames (%.1f s)',...
              idx, fps, percentile_window, percentile_window/fps));
grid on;
set(ax1, 'TickLength', [0 0]);

ax2 = subplot(312);
plot(t, tr2, 'Color', c2);
hold on;
plot(t, b2, 'k-');
hold off;
ylabel('Least squares');
grid on;
set(ax2, 'TickLength', [0 0]);

ax3 = subplot(313);
plot(t, tr1_dff, 'Color', c1);
hold on;
plot(t, tr2_dff, '--', 'Color', c2);
hold off;
ylabel('DFF');
grid on;
set(ax3, 'TickLength', [0 0]);

linkaxes([ax1 ax2 ax3], 'x');
xlim(t([1 end]));
xlabel('Time (s)');

end

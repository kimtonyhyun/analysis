function [tr2_dff, info2, ax2] = inspect_dff_traces(ds_proj, ds_ls, idx, fps, varargin)
% DFF traces are most straightforward to compute when using simple
% projection to the motion-corrected imaging data (with minimal additional
% processing). However, least-squares projection provides better safeguards
% against trace contamination.
%
% This basic visualization function compares DFF traces obtained by simple
% projection vs. least-squares, as a sanity check that the DFF scaling of
% the two methods are comparable. In addition, the baseline (F0) estimate
% is shown for both projection and least squares traces. The 'varargin'
% parameters are passed to 'compute_dff_traces', which in turn passes those
% to 'fix_baseline'.
%
% Returns: Results of 'compute_dff_trace' for the least squares trace
%

% Default colors
c1 = [0 0.447 0.741];
c2 = [0.466 0.674 0.188];

% Custom "subplot" command that leaves less unusued space between panels
sp = @(m,n,p) subtightplot(m, n, p, 0.05, 0.05, 0.05); % Gap, Margin-X, Margin-Y

tr1 = ds_proj.get_trace(idx);
[tr1_dff, info1] = compute_dff_trace(tr1, varargin{:});
b1 = info1.baseline;

tr2 = ds_ls.get_trace(idx);
[tr2_dff, info2] = compute_dff_trace(tr2, varargin{:});
b2 = info2.baseline;
nu2 = calculate_noise_level_nu(tr2_dff, fps);

num_frames = length(tr1);
t = 1/fps*(0:num_frames-1);

ax1 = sp(3,1,1);
switch info1.method.name
    case 'nonactive_polyfit'
        draw_nonactive_polyfit(t, tr1, info1, c1);
        
    otherwise
        plot(t, tr1, 'Color', c1);
end
hold on;
plot(t, b1, 'k-', 'LineWidth', 2);
hold off;
ylabel('Simple projection');
title(sprintf('Cell %d: FPS=%.1f Hz', idx, fps));
grid on;
ylim(compute_y_range(tr1));
set(ax1, 'TickLength', [0 0]);

ax2 = sp(3,1,2);
switch info2.method.name
    case 'nonactive_polyfit'
        draw_nonactive_polyfit(t, tr2, info2, c2);
        
    otherwise
        plot(t, tr2, 'Color', c2);
end
hold on;
plot(t, b2, 'k-', 'LineWidth', 2);
hold off;
ylabel('Least squares');
grid on;
ylim(compute_y_range(tr2));
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

end % inspect_dff_traces

function y_range = compute_y_range(tr)
    m = min(tr);
    M = max(tr);
    y_range = [m M] + 0.1*(M-m)*[-1 1];
end

function draw_nonactive_polyfit(t, tr, info, color)

cla; hold on;

nf = info.method.nonactive_frames;
af = ~nf; % active_frames

nf_segs = frame_list_to_segments(find(nf));
for k = 1:size(nf_segs,1)
    frames = nf_segs(k,1):nf_segs(k,2);
    plot(t(frames), tr(frames), 'Color', color);
end

af_segs = frame_list_to_segments(find(af));
for k = 1:size(af_segs,1)
    frames = af_segs(k,1):af_segs(k,2);
    plot(t(frames), tr(frames), 'r');
end

tr_thresh = info.method.tr_thresh;
plot(t([1 end]), tr_thresh*[1 1], 'k--');
hold off;

end % draw_nonactive_polyfit

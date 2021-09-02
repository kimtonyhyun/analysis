function [trace_out, info] = fix_baseline(trace_in, varargin)
% Tracks wandering baseline of a fluorescence trace by computing the
% running percentile. Previously used in Wagner et al. Cell, 2019.
%
% Be aware that, for highly active Ca2+ traces, the running percentile may
% not be the best way of estimating the "baseline"!

% Default parameters. As used in Wagner et al. 2019.
p = 10; % 10-th percentile
window = 1000; % In frames. At 30 fps, this is ~33.3 s
compute_dff = false;

for k = 1:length(varargin)
    if ischar(varargin{k})
        switch lower(varargin{k})
            case {'p', 'percentile'}
                p = varargin{k+1};
            case 'window'
                window = varargin{k+1};
            case 'dff'
                compute_dff = true;
        end
    end
end

trace_shape = size(trace_in);
trace_in = trace_in(:);

baseline = running_percentile(trace_in, window, p);
baseline = correct_bias(trace_in, baseline);
trace_out = trace_in - baseline;
if compute_dff
    if any(baseline < 0)
        cprintf('red', 'Warning: F0 trace contains negative values\n');
    end
    trace_out = (trace_in - baseline)./baseline;
end
trace_out = reshape(trace_out, trace_shape);

info.name = 'fix_baseline';
info.percentile = p;
info.window = window;
info.baseline = reshape(baseline, trace_shape);

end % fix_baseline

function baseline2 = correct_bias(trace, baseline)
% Baseline estimation techniques that fits to the lower percentile data
% points of a trace typically has a negative bias. This function attempts
% to correct that bias.
    stats = compute_trace_stats(trace - baseline);
    baseline2 = baseline + stats.mode;
end


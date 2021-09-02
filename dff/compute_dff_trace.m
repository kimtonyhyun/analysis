function [tr_dff, info] = compute_dff_trace(tr, varargin)
% Minor wrapper around 'fix_baseline'. (TODO: Consider merging?)

[df_trace, info] = fix_baseline(tr, varargin{:});

if any(info.baseline < 0)
    cprintf('blue', '  Warning: Baseline (F0) is negative\n');
end

tr_dff = df_trace./info.baseline;
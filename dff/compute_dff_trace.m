function [tr_dff, info] = compute_dff_trace(tr, varargin)
% To compute "proper" DFF traces, the fluorescence trace must be scaled
% such that a value of "0" indicates the optical signal in the absence of 
% fluorescence (i.e. measurement with the excitation light off).

% Estimate the baseline for DFF computation
[~, fix_info] = fix_baseline(tr, varargin{:});
baseline = fix_info.baseline;

if any(baseline < 0)
    cprintf('blue', '  Warning: Baseline (F0) is negative\n');
end

% Prep output
tr_dff = (tr - baseline)./baseline;

info.fix_params = fix_info.method;
info.baseline = baseline;
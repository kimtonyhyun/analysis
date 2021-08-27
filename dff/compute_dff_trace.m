function [tr_dff, info] = compute_dff_trace(tr, varargin)
% To compute "proper" DFF traces, the fluorescence trace must be scaled
% such that a value of "0" indicates the optical signal in the absence of 
% fluorescence (i.e. measurement with the excitation light off).

% Defaults
fix_wandering_baseline = false;
percentile_window = 1e3; % Frames

for k = 1:length(varargin)
    if ischar(varargin{k})
        switch lower(varargin{k})
            case {'fix', 'fix_baseline'}
                fix_wandering_baseline = true;
                percentile_window = varargin{k+1};
        end
    end
end

if fix_wandering_baseline
    % Estimate the baseline, which may change slowly over the recording.
    % Track the wandering baseline by calculating the 10th percentile trace
    % value over the specified number of frames.
    [~, fix_info] = fix_baseline(tr, 'percentile', 'window', percentile_window);
    baseline = fix_info.baseline;

    % Correct vertical offset originating from the use of the 10th percentile
    stats = compute_trace_stats(tr - baseline);
    offset = stats.mode;
    baseline = baseline + offset; 
    
    info.fix_params = fix_info.method;
else
    % Baseline is estimated as a scalar
    stats = compute_trace_stats(tr);
    F0 = stats.mode;
    baseline = F0 * ones(size(tr), 'single');
end

% Prep output
tr_dff = (tr - baseline)./baseline;

info.fix_wandering_baseline = fix_wandering_baseline;
info.baseline = baseline;
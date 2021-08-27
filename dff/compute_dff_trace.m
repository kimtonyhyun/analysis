function [tr_dff, info] = compute_dff_trace(tr, varargin)
% To compute "proper" DFF traces, the fluorescence trace must be scaled
% such that a value of "0" indicates the optical signal in the absence of 
% fluorescence (i.e. measurement with the excitation light off).

% Defaults
percentile_window = 3e3; % Frames

for k = 1:length(varargin)
    if ischar(varargin{k})
        switch lower(varargin{k})
            case 'window'
                percentile_window = varargin{k+1};
        end
    end
end

% Estimate the baseline (F0), which may change slowly over the recording.
% We track the wandering baseline by calculating the 10th percentile trace
% value over the specified number of frames.
[~, fix_info] = fix_baseline(tr, 'percentile', 'window', percentile_window);
baseline = fix_info.baseline;

% Correct vertical offset originating from the use of the 10th percentile
stats = compute_trace_stats(tr - baseline);
offset = stats.mode;
baseline = baseline + offset;

% Prep output
tr_dff = (tr - baseline)./baseline;

info.fix_params = fix_info.method;
info.baseline = baseline;
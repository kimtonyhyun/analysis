function [trace_out, info] = fix_baseline(trace_in, method_name, varargin)

% Default method (to be consistent with past behavior)
if ~exist('method_name', 'var')
    method_name = 'percentile';
end

trace_shape = size(trace_in);
trace_in = trace_in(:);

switch (method_name)
    % Each method must set 'trace_out', 'baseline', and 'method'
    %------------------------------------------------------------
    case 'mode'
        stats = compute_trace_stats(trace_in);
        F0 = stats.mode;
        baseline = F0 * ones(size(trace_in), 'single');
        trace_out = trace_in - baseline;
        
        method.name = 'mode';
    
    case {'linear', 'lin'}
        % Default parameters
        % TODO: Implement
        
    case {'percentile', 'prctile'}
        % Default parameters. As used in Wagner et al. Cell (2019).
        p = 10; % 10-th percentile
        window = 1000; % In frames. At 30 fps, this is ~33.3 s
        
        for k = 1:length(varargin)
            if ischar(varargin{k})
                switch lower(varargin{k})
                    case 'p'
                        p = varargin{k+1};
                    case 'window'
                        window = varargin{k+1};
                end
            end
        end
        
        baseline = running_percentile(trace_in, window, p);
        baseline = correct_offset(trace_in, baseline);
        trace_out = trace_in - baseline;
        
        method.name = 'percentile';
        method.percentile = p;
        method.window = window;
end
 
info.method = method;
info.baseline = reshape(baseline, trace_shape);

trace_out = reshape(trace_out, trace_shape);

end % fix_baseline

function baseline2 = correct_offset(trace, baseline)
% Baseline estimation techniques that fits to the lower percentile data
% points of a trace typically has a bias to the "true" baseline. This
% function attempts to correct that negative bias.
    stats = compute_trace_stats(trace - baseline);
    baseline2 = baseline + stats.mode;
end
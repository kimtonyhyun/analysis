function [trace_out, info] = fix_baseline(trace_in, method_name, varargin)
% Note: Baseline is subtracted off the trace. Thus, one should _not_ run
%   fix_baseline when accurate DFF values are desired.
%
if ~exist('method_name', 'var')
    method_name = 'percentile';
end

trace_shape = size(trace_in);
trace_in = trace_in(:);

switch (method_name)       
    case {'percentile', 'prctile'}
        % Default parameters
        p = 10; % 10-th percentile
        moving_window_num_frames = 1000; % At 30 fps, this is ~33.3 s
        
        for k = 1:length(varargin)
            if ischar(varargin{k})
                switch lower(varargin{k})
                    case 'p'
                        p = varargin{k+1};
                    case 'window'
                        moving_window_num_frames = varargin{k+1};
                end
            end
        end
        
        baseline = running_percentile(trace_in, moving_window_num_frames, p);
        trace_out = trace_in - baseline;
        
        method.name = 'percentile';
        method.p = p;
        method.moving_window_num_frames = moving_window_num_frames;
end
 
info.method = method;
info.baseline = reshape(baseline, trace_shape);

end % fix_baseline
function [trace_out, info] = fix_baseline(trace_in, method_name, varargin)

if ~exist('method_name', 'var')
    method_name = 'hakan'; % Legacy
end

trace_in = trace_in(:);

switch (method_name)
    case 'hakan'
        trace_out = fix_baseline_hakan(trace_in);
        
        method.name = 'hakan';
        
    case {'percentile', 'prctile'}
        % Default parameters
        p = 10; % 10-th percentile
        fps = 30; % Frames per second
        moving_window_duration = 30; % s
        
        for k = 1:length(varargin)
            if ischar(varargin{k})
                switch lower(varargin{k})
                    case 'p'
                        p = varargin{k+1};
                    case 'fps'
                        fps = varargin{k+1};
                    case 'window'
                        moving_window_duration = varargin{k+1};
                end
            end
        end
        
        baseline = running_percentile(trace_in, fps * moving_window_duration, p);
        trace_out = trace_in - baseline;
        
        method.name = 'percentile';
        method.p = p;
        method.fps = fps;
        method.moving_window_duration = moving_window_duration;
end
 
info.method = method;

end % fix_baseline
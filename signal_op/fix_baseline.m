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
    
    case 'nonactive_polyfit'
        % Default parameters
        
        % For identifying non-active portions of the trace
        thresh = 0.5;
        padding = 100; % In frames
        
        % For polyfit
        order = 1;
        
        for k = 1:length(varargin)
            if ischar(varargin{k})
                switch lower(varargin{k})
                    case 'thresh'
                        thresh = varargin{k+1};
                    case 'padding'
                        padding = varargin{k+1};
                    case 'order'
                        order = varargin{k+1};
                end
            end
        end
        
        [baseline, tr_thresh, nonactive_frames] = fit_nonactive_frames(trace_in, thresh, padding, order);
        trace_out = trace_in - baseline;
        
        method.name = 'nonactive_polyfit';
        method.thresh = thresh;
        method.padding = padding;
        method.order = order;
        method.tr_thresh = tr_thresh;
        method.nonactive_frames = nonactive_frames;
        
    case 'percentile'
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

function [baseline, tr_thresh, nonactive_frames] = fit_nonactive_frames(trace, thresh, padding, order)

    num_frames = length(trace);
    t = linspace(0, 1, num_frames)';

    tr_min = min(trace);
    tr_max = max(trace);
    tr_thresh = tr_min + thresh * (tr_max - tr_min);

    active_frames = parse_active_frames(trace > tr_thresh, padding);
    active_frames = frame_segments_to_list(active_frames);
    
    nonactive_frames = true(num_frames, 1);
    nonactive_frames(active_frames) = false;

    if (sum(nonactive_frames) < order)
        cprintf('red',  '  Insufficient number of nonactive frames for baseline fitting!\n');
        tr_thresh = tr_max;
        nonactive_frames = true(num_frames, 1);
    end
    
    x = t(nonactive_frames);
    y = trace(nonactive_frames);

    coeffs = polyfit(x, y, order);
    baseline = polyval(coeffs, t);
end
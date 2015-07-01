function idx_events = extract_events(trace,varargin)

% Extract local maxima in the trace as events
%
% trace: vector of fluroescence values to event-detect from
%
% Variable input arguments:
%
% 'visualize' : add it as an argument to plot the output
% 'mad_scale' : enter a numeric value following it to set custom scale for
%               the mad (median absolute deviation) of the trace for 
%               determining an event threshold
%
% Output:
% idx_events: indices containing the detected events
%

% Defaults
mad_scale = 8; 
do_plot = 0;

if ~isempty(varargin)
    for k = 1:length(varargin)
        switch varargin{k}
            case 'mad_scale'
                mad_scale = varargin{k+1};
                if ~isnumeric(mad_scale)
                    error('Mad scale must be numerical');
                end
            case 'visualize'
                do_plot = 1;
        end
    end
end

num_frames = length(trace);
smooth_trace = medfilt1(trace,5);
mad = compute_mad(smooth_trace);
thresh = mad_scale * mad;

idx_events = calc_localmax_above_threshold(thresh,smooth_trace);
val_events = zeros(1,length(idx_events));

%go back to the original trace and spot the exact maxima
for k = 1:length(idx_events)
    coarse_event = idx_events(k);
    go_back = min(5,coarse_event-1);
    go_forward = min(5,num_frames-coarse_event);
    [maxval,idx_max] = max(trace((coarse_event-go_back):(coarse_event+go_forward)));
    idx_events(k) = idx_max+coarse_event-go_back-1;
    val_events(k) = maxval;
end

% Remove duplicates (may happen because of median filtering)
[idx_events,idx_unique] = unique(idx_events);
val_events = val_events(idx_unique);

if do_plot    
    spikes_vec = zeros(1,length(smooth_trace));        
    spikes_vec(idx_events) = val_events;

    plot(smooth_trace);
    hold on;
    stem(spikes_vec);
    plot(ones(1,length(smooth_trace))*thresh,'r')
    plot(trace,'--m');
    hold off
    for k = 1:length(idx_events)
        text(idx_events(k),double(thresh*1.1),num2str(k))
    end
end
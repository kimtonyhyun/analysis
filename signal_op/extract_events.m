function event_times = extract_events(trace,varargin)
% function event_times = extract_events_2step(cell_of_interest)
%
% Extract event times from the ICA trace of a cell doing 2 iterations
%
%   cell_of_interest : Nx1 cell array where each entry contains the ICA trace for
%   one trial, with total of N trials
% 
%   event_times : Nx1 cell where each entry is a vector whose elements are
%   the frame indices where an event occurs
%
% Hakan Inan (Nov 2014)
%

%%%%
%some parameters that can be changed within the script (advanced):
%%%%
%mad for df/f trace is multiplied with below to get the coarse threshold
coarse_thresh_multiplier = 9;
%mad for deconvolved trace is multiplied with below to get the fine threshold
fine_thresh_multiplier = 15;
%offset of the impulse for deconvolution
offset = 0;

% Length before and after the maximum of the characteristic event, h
num_frames_before_max = 20;
num_frames_after_max = 100;


idx_event_frames = extract_candidate_events(trace,coarse_thresh_multiplier,0);
out = iterate_event_indices(idx_event_frames);
idx_event_frames_v2 = out{1};h_1stround = out{2};
out = iterate_event_indices(idx_event_frames_v2);
idx_event_frames_v3 = out{1};h_2ndround = out{2};g_2ndround = out{3};
event_times = idx_event_frames_v3;


%%%%%%%%%%%%%%
%visualization
%%%%%%%%%%%%%%
if ~isempty(varargin)
    if varargin{1}==1
        figure,
        plot(h_1stround,'k');
        hold on
        plot(h_2ndround,'r');
        hold off
        legend('the characteristic event, 1st iteration','the characteristic event, 2nd iteration');
        xlabel('frame')

        trace = trace / max(trace);
        mad_scale_coarse = compute_mad(trace);
        threshold_coarse = coarse_thresh_multiplier*mad_scale_coarse;

        filtered =  conv(trace,g_2ndround,'same');
        filtered = filtered / max(filtered);
        mad_scale_fine = compute_mad(filtered);
        threshold_fine = fine_thresh_multiplier*mad_scale_fine;

        spikes_vec_fine = zeros(1,length(trace));
        spikes_vec_fine(idx_event_frames_v3) = 1;
        spikes_vec_coarse = zeros(1,length(trace));
        spikes_vec_coarse(idx_event_frames) = 1;

        figure,
        plot(trace,'k');
        hold on
        plot(filtered,'m');
        plot(1:length(trace),repmat(threshold_coarse,1,length(trace)),'--b');%coarse threshold
        plot(1:length(trace),repmat(threshold_fine,1,length(trace)),'--r');%fine threshold
        stem(spikes_vec_fine*threshold_fine);
        stem(spikes_vec_coarse*threshold_coarse);
        hold off
        xlabel('frame')
        ylabel('a.u.')
        title('raw and deconvolved ica traces together with the thresholds and the detected events')
        legend('raw trace','deconvolved trace','coarse threshold','deconv threshold')
    end
end

%%%%%%%%%%%%%%%%%%%
%internal functions
%%%%%%%%%%%%%%%%%%%
function out = iterate_event_indices(idx_event_frames_in)

    h_len = num_frames_before_max+num_frames_after_max+1;
    %surrogate waveform
    event_mat = construct_events(idx_event_frames_in);
    %calculate the mean of all the trials in the event_mat discarding zeros
    h = zeros(h_len,1);
    for k = 1:h_len
        dum = event_mat(:,k);
        nonzeros = sum(dum~=0);
        if nonzeros>=1
            h(k) = sum(dum)/nonzeros;
        end
    end
    %optional: smooth h before inverting
    SGfilter = SGsmoothingfilter(20,8);
    h = conv(h,SGfilter,'same');
    g_len= h_len; %try to find an inverse g that has the same length as h
    g = invert_h(h,g_len); %invert h 
    %convolve g with the trace of the whole day and use it to calculate an
    %adaptive threshold for detecting events
    trace_filtered = conv(trace,g,'same');
    idx_event_frames_out = extract_candidate_events(trace_filtered,fine_thresh_multiplier,0);
    out = {idx_event_frames_out,h,g};
end

function event_mat = construct_events(idx_event_frames)

%align all the trial segments(where there is an event) and accumulate all of them.
%The events in event_mat are of fixed length with frames before the maximum and
%frames after the maximum determind by the function inputs.
%If an event is not long enough on either side, use all the frames to the
%end on the respective side and append zeros.

    h_len = num_frames_before_max+num_frames_after_max+1;
    event_mat = zeros(length(idx_event_frames),h_len);

    for k = 1:length(idx_event_frames) 

        trace_len = length(trace);
        idx_max = idx_event_frames(k);

        if k<length(idx_event_frames) % k is not the last event
            num_frames_to_next = idx_event_frames(k+1)-idx_max-1;
            numframes2add_after = max(min(num_frames_after_max,num_frames_to_next-10),1);
        else % k is the last event
            numframes2add_after = min(trace_len-idx_max,num_frames_after_max);
        end
        if k>1 % k is not the first event
            num_frames_from_previous = idx_max-idx_event_frames(k-1);
            numframes2add_before = max(min(num_frames_before_max,num_frames_from_previous-10),1);
        else % k is the first event
            numframes2add_before = min(idx_max-1,num_frames_before_max);
        end

        indices_trace = (-numframes2add_before:numframes2add_after)+idx_max;
        indices_event = (-numframes2add_before:numframes2add_after) + num_frames_before_max+1;

        event_to_add = trace(indices_trace);
        event_mat(k,indices_event) = event_to_add/max(event_to_add);
    end

end

function g = invert_h(h,g_len)
    h_len = numel(h);
    Mconv =  toeplitz([h;zeros(g_len-1,1)],zeros(1,g_len));
    t_desired = floor(g_len/2)+num_frames_before_max+offset;
    y = zeros(h_len+g_len-1,1);
    y(t_desired) = 1;

    % Append noise to A for stability
    noise_length = ceil(h_len)*4;
    noise = wgn(noise_length,g_len,0)*0.4; 

    A = [Mconv;noise];
    y = [y;zeros(noise_length,1)];
    
    % Below commented out portion does a better job at filter inversion,
    % however runs slower and requires CVX package
    % A_not = A([1:(t_desired-1),(t_desired+1):end],:);
    % cvx_begin quiet
    %     variable g(g_len,1)
    %     minimize(norm(A_not*g,1))
    %     subject to
    %         A(t_desired,:)*g == 1;
    % cvx_end

    g = pinv(A)*y;
end

end
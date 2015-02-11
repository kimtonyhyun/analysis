function event_mat = construct_events_from_max(cell_of_interest,...
    max_mat,frames_before_max,frames_after_max)
%align all the trial segments(where there is an event) and accumulate all of them.
%The events in event_mat are of fixed length with frames before the maximum and
%frames after the maximum determind by the function inputs.
%If a trial is not long enough on either side, use all the frames to the
%end on the respective side and append zeros.
num_trials_with_events = length(unique(max_mat(:,1)));
h_len = frames_before_max+frames_after_max+1;
event_mat = zeros(num_trials_with_events,h_len);
cnt=0; %counter for the loop 
trials_with_events = max_mat(:,1)';
for k = trials_with_events 
    cnt = cnt+1;
    trial = cell_of_interest{k};
    %handle overlapping waveforms within a single trial
    if cnt < num_trials_with_events
        if trials_with_events(cnt+1) == k %overlapping waveform
            idx_maxi_this = max_mat(cnt,2);
            idx_maxi_that = max_mat(cnt+1,2);
            if idx_maxi_that-idx_maxi_this < frames_before_max+frames_after_max
                %cut the trial before the other event begins
                trial((idx_maxi_that-frames_before_max):end) = [];
            end
        end
    end
    trial_length = length(trial);
    idx_maxi = max_mat(cnt,2);
    maxval = max(trial);
    if idx_maxi<frames_before_max+1 %we have less frames to the left
        idx_start_trial = 1;
        idx_start_event = 1+ (1+frames_before_max-idx_maxi);
    else 
        idx_start_trial = idx_maxi-frames_before_max;
        idx_start_event = 1;
    end
    if trial_length-idx_maxi<frames_after_max %we have less frames to the right
        idx_end_trial = trial_length;
        idx_end_event = h_len - (frames_after_max - (trial_length-idx_maxi) );
    else
        idx_end_trial = idx_maxi+frames_after_max;
        idx_end_event = h_len;
    end
    event_mat(cnt,idx_start_event:idx_end_event) = trial(idx_start_trial:idx_end_trial)/maxval;
end
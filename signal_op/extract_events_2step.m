function event_times = extract_events_2step(cell_of_interest)
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
%some parameters that can be changed within the script:
%%%%
%std for df/f trace is multiplied with below to get the coarse threshold
coarse_thresh_multiplier = 9;
%std for deconvolved trace is multiplied with below to get the fine threshold
fine_thresh_multiplier = 15;
%offset of the impulse for deconvolution
offset = 40;

num_trials = length(cell_of_interest);
cell_concat = [];

%construct a cell of traces for all the trials of the given session
for k = 1:num_trials
    cell_concat = [cell_concat,cell_of_interest{k}];
end
%compute a coarse threshold for detecting event in any trial:
mad = compute_mad(cell_concat);
thresh_coarse = mad*coarse_thresh_multiplier;
%compute the trials with events:
max_mat = [];
for k = 1:num_trials
    [maxi,idx_maxi] = max(cell_of_interest{k});
    if maxi>thresh_coarse
        max_mat = [max_mat;[k,idx_maxi]];
    end
end
%construct a matrix of events to form a 'surrogate' waveform from 
frames_before_max = 20;
frames_after_max = 100;
h_len = frames_before_max+frames_after_max+1;
event_mat = construct_events_from_max(cell_of_interest,...
    max_mat,frames_before_max,frames_after_max);
%plot(event_mat');
%calculate the mean of all the trials in the event_mat discarding zero
%terms
h = zeros(1,h_len);
for k = 1:h_len
    dum = event_mat(:,k);
    nonzeros = sum(dum~=0);
    if nonzeros>1
        h(k) = sum(dum)/nonzeros;
    end
end
h_1stround = h;
%optional: smooth h before inverting
h = SG_filter_trace(h,20,8);
g_len= h_len; %try to find an inverse g that has the same length as h
g = invert_filter(h,g_len,frames_before_max+offset);

%convolve g with the trace of the whole day and use it to calculate an
%adaptive threshold for detecting events
cell_filtered = conv(cell_concat,g);
mad = compute_mad(cell_filtered);
threshold = mad*fine_thresh_multiplier; 

%convolve g with the traces to obtain event times for each trial
event_times = cell(num_trials,1);
for k = 1:num_trials
    out = conv(cell_of_interest{k},g);
    out = out(1+offset:end-offset);
    out_smooth =  SG_filter_trace(out,10,8);
    event_times{k} = calc_localmax_above_threshold(threshold,out_smooth);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% second round of iteration:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_mat2 = [];
for k = 1:num_trials
    event_time_vec = event_times{k};
    if ~isempty(event_time_vec)
        for l = 1:length(event_time_vec)
            max_mat2 = [max_mat2;[k,event_time_vec(l)]];
        end
    end
end
%it might be that there are no events in this iteration. In that case
%simply use the max_mat from the previous iteration
if ~isempty(max_mat2)
    max_mat = max_mat2;
end    

%surrogate waveform
h_len = frames_before_max+frames_after_max+1;
event_mat = construct_events_from_max(cell_of_interest,...
    max_mat,frames_before_max,frames_after_max);
%calculate the mean of all the trials in the event_mat discarding zeros
h = zeros(1,h_len);
for k = 1:h_len
    dum = event_mat(:,k);
    nonzeros = sum(dum~=0);
    if nonzeros>=1
        h(k) = sum(dum)/nonzeros;
    end
end
h_2ndround = h;
%optional: smooth h before inverting
h = SG_filter_trace(h,20,8);
g_len= h_len; %try to find an inverse g that has the same length as h
g = invert_filter(h,g_len,frames_before_max+offset); %invert h 
%convolve g with the trace of the whole day and use it to calculate an
%adaptive threshold for detecting events
cell_filtered = conv(cell_concat,g);
mad = compute_mad(cell_filtered);
threshold = mad*fine_thresh_multiplier; 
%convolve g with the traces to obtain event times for each trial
event_times = cell(num_trials,1);
for k = 1:num_trials
    out = conv(cell_of_interest{k},g);
    out = out(1+offset:end-offset);
    out_smooth =  SG_filter_trace(out,10,8);
    event_times{k} = calc_localmax_above_threshold(threshold,out_smooth);
end

%%%%%%%%%%%%%%
%visualization
%%%%%%%%%%%%%%
figure,
plot(h_1stround,'k');
hold on
plot(h_2ndround,'r');
hold off
legend('the surrogate event, 1st iteration','the surrogate event, 2nd iteration');
xlabel('frame')
cell_concat = SG_filter_trace(cell_concat,10,8);
normalizer_coarse = max(cell_concat);
cell_concat = cell_concat / normalizer_coarse;
thresh_coarse = thresh_coarse / normalizer_coarse;
event_times_coarse = calc_localmax_above_threshold(thresh_coarse,cell_concat);

cell_filtered = cell_filtered(1+offset:end-offset);
out_smooth = SG_filter_trace(cell_filtered,10,8);
normalizer_fine = max(out_smooth);
out_smooth = out_smooth / normalizer_fine;
thresh_fine = threshold / normalizer_fine;
event_times_fine = calc_localmax_above_threshold(thresh_fine,out_smooth);

spikes_vec_fine = zeros(1,length(cell_concat));
spikes_vec_fine(event_times_fine) = 1;
spikes_vec_coarse = zeros(1,length(cell_concat));
spikes_vec_coarse(event_times_coarse) = 1;

figure,
plot(cell_concat,'k');
hold on
plot(out_smooth,'m');
plot(1:length(cell_concat),repmat(thresh_coarse,1,length(cell_concat)),'--b');%coarse threshold
plot(1:length(cell_concat),repmat(thresh_fine,1,length(cell_concat)),'--r');%fine threshold
stem(spikes_vec_fine*thresh_fine);
stem(spikes_vec_coarse*thresh_coarse);
hold off
xlabel('frame')
ylabel('a.u.')
title('smoothed and deconvolved ica traces together with the thresholds and the detected events')
legend('smoothed trace','deconvolved trace','coarse threshold','deconv threshold')
%}
% out = [];
% for k = 1:size(event_mat,1)
%     in = event_mat(k,:);
%     inp = SG_filter_trace(in,20,8);
%     dum = conv(g,inp);
%     dum = dum(floor(g_len/2)+1:end-floor(g_len/2));
%     out = [out;dum];
% end
% plot(out')

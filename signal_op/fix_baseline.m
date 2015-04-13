function trace_out = fix_baseline(trace_in)
% Remove baseline offset of cell trace by fitting a line to it and
% subtracting

num_frames = length(trace_in);
h = polyfit(1:num_frames,min(trace_in,2*quantile(trace_in,0.3)),1);
subst = h(2)+h(1)*(1:num_frames);
trace_out = trace_in-subst;
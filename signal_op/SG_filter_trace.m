function out = SG_filter_trace(in,N,L)
h = SGsmoothingfilter(N,L);
out = conv(h,in); %apply SG filter on the raw data
out = out(N+1:end-N); %leave out the transients of the filtering



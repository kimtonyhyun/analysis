function out = SG_filter_trace(in,N,L)
 h = SGsmoothingfilter(N,L);
out = conv(in,h,'same'); %apply SG filter on the raw data




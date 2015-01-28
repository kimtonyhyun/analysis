function save_classification(class, outfile)

num_traces = length(class);

fid = fopen(outfile, 'w');
for trace_idx = 1:num_traces
    fprintf(fid, '%d, %s\n', trace_idx, class{trace_idx});
end
fclose(fid);
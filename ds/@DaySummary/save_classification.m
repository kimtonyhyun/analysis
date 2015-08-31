function save_classification(obj, outfile)

num_candidates = obj.num_cells;

fid = fopen(outfile, 'w');
for k = 1:num_candidates
    fprintf(fid, '%d, %s\n', k, obj.cells(k).label);
end
fclose(fid);
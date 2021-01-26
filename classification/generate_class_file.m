function class_savename = generate_class_file(num_cells)
% Generate a "dummy" classification file, where all 'num_cells' cells are
% specified to be a cell

timestamp = datestr(now, 'yymmdd-HHMMSS');
class_savename = sprintf('class_%s.txt', timestamp);

fid = fopen(class_savename, 'w');
for k = 1:num_cells
    fprintf(fid, '%d, cell\n', k);
end
fclose(fid);
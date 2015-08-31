function load_class(obj, source)

fid = fopen(source, 'r');
class = textscan(fid, '%d %s', 'Delimiter', ',');
fclose(fid);

class = class{2};
for k = 1:obj.num_cells
    obj.cells(k).label = class{k};
end
function class = load_classification(source)

fid = fopen(source, 'r');
class = textscan(fid, '%d %s', 'Delimiter', ',');
fclose(fid);

class = class{2};
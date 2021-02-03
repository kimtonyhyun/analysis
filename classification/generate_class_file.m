function class_savename = generate_class_file(num_cells, varargin)
% Generate a "dummy" classification file, where all 'num_cells' cells are
% specified to be a cell

timestamp = datestr(now, 'yymmdd-HHMMSS');

for k = 1:length(varargin)
    vararg = varargin{k};
    if ischar(vararg)
        switch lower(vararg)
            case 'timestamp'
                timestamp = varargin{k+1};
        end
    end
end

class_savename = sprintf('class_%s.txt', timestamp);

fid = fopen(class_savename, 'w');
for k = 1:num_cells
    fprintf(fid, '%d, cell\n', k);
end
fclose(fid);
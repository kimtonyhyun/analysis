function [pos, neg] = collect_tdt_labels(ds)

num_cells = ds.num_cells;

pos = zeros(1, num_cells); % Preallocate
neg = zeros(1, num_cells);

num_pos = 0;
num_neg = 0;

for k = 1:num_cells
    tdt_label = ds.cells(k).label;
    if ~isempty(tdt_label)
        switch tdt_label
            case 'positive'
                num_pos = num_pos + 1;
                pos(num_pos) = k;

            case 'negative'
                num_neg = num_neg + 1;
                neg(num_neg) = k;
        end
    end
end

pos = pos(1:num_pos);
neg = neg(1:num_neg);
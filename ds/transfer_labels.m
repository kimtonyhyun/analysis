function [transfer_to_1, transfer_to_2] = transfer_labels(ds1, ds2, m_1to2, m_2to1)

% Transfer (positive) labels from ds2 --> ds1
%------------------------------------------------------------
is_cell1 = find(ds1.is_cell);
transfer_to_1 = compute_ids_to_transfer(ds1, ds2, m_1to2);
num_transfer_to_1 = length(transfer_to_1);
is_not_cell1 = setdiff(find(~ds1.is_cell), transfer_to_1);

subplot(121);
ds1.plot_cell_map({is_cell1, 'g'; transfer_to_1, 'c'; is_not_cell1, 'r'});
title(sprintf('Dataset 1: %d transferred labels', num_transfer_to_1));

% Transfer (positive) labels from ds1 --> ds2
%------------------------------------------------------------
is_cell2 = find(ds2.is_cell);
transfer_to_2 = compute_ids_to_transfer(ds2, ds1, m_2to1);
num_transfer_to_2 = length(transfer_to_2);
is_not_cell2 = setdiff(find(~ds2.is_cell), transfer_to_2);

subplot(122);
ds2.plot_cell_map({is_cell2, 'g'; transfer_to_2, 'c'; is_not_cell2, 'r'});
title(sprintf('Dataset 2: %d transferred labels', num_transfer_to_2));

% Update the labels in place
%------------------------------------------------------------
fprintf('  Press any key to update labels in place\n');
pause;
for cell_idx1 = transfer_to_1
    ds1.cells(cell_idx1).label = 'cell';
end
for cell_idx2 = transfer_to_2
    ds2.cells(cell_idx2).label = 'cell';
end

end % transfer_labels

function transfer_ids = compute_ids_to_transfer(ds_target, ds_source, m_ttos)
    % Find cells in 'ds_target' that are not cells, but their corresponding
    % cell in 'ds_source' are classified cells.
    is_not_cell = find(~ds_target.is_cell); % Indices of non-cells in target
    is_cell_in_source = false(size(is_not_cell));

    for k = 1:length(is_not_cell)
        cell_idx1 = is_not_cell(k);
        match = m_ttos{cell_idx1};
        if ~isempty(match)
            cell_idx2 = match(1);
            is_cell_in_source(k) = ds_source.is_cell(cell_idx2);
        end
    end
    transfer_ids = is_not_cell(is_cell_in_source);
end % transfer_labels
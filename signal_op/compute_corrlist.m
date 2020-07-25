function corrlist = compute_corrlist(ds1, ds2)
% Assumes that the traces of the two DaySummary's have been resampled so
% that they have the same time base.
%
% Note that it is possible to input the same DaySummary 

cell_inds1 = find(ds1.is_cell);
traces1 = ds1.get_trace(cell_inds1)'; % [num_samples x num_cells]

if ~exist('ds2', 'var') % Only ds1 provided
    corrlist = corr_to_corrlist(corr(traces1), 'upper');
    
    corrlist(:,1) = cell_inds1(corrlist(:,1));
    corrlist(:,2) = cell_inds1(corrlist(:,2));
else % ds2 also provided
    cell_inds2 = find(ds2.is_cell);
    traces2 = ds2.get_trace(cell_inds2)';
    corrlist = corr_to_corrlist(corr(traces1, traces2));
    
    corrlist(:,1) = cell_inds1(corrlist(:,1));
    corrlist(:,2) = cell_inds2(corrlist(:,2));
end
corrlist = sortrows(corrlist, 3, 'descend');

function [boundary, ic_mask] = compute_ic_boundary(ic_filter, ic_filter_thresh)
% Computes the boundary of the IC filter

A = threshold_ic_filter(ic_filter, ic_filter_thresh);

ic_mask = (A>0);
boundaries = bwboundaries(ic_mask, 'noholes');
boundary = boundaries{1};
for i = 2:length(boundaries)
    new_boundary = boundaries{i};
    if (size(new_boundary,1) > size(boundary,1))
        boundary = new_boundary;
    end
end

boundary = fliplr(boundary); % Result arranged as [x y]
function filter_out = rescale_filter_to_clim(filter, clim)
% Numerically rescale the filter to match the provided clim

    f_max = max(filter(:));
    f_min = min(filter(:));

    filter_norm = (filter-f_min)/f_max; % Matched to range [0, 1]

    clim_delta = clim(2)-clim(1);
    clim_usage = 0.8;
    filter_out = clim(1) + (1-clim_usage)/2*clim_delta +...
                 clim_usage*clim_delta*filter_norm;

end
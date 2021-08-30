function nu = calculate_noise_level_nu(trace, fps)
% Calculate the DF/F noise level per Rupprecht et al. 2021

abs_diff_tr = abs(diff(trace));
nu = median(abs_diff_tr) / sqrt(fps);
nu = 100 * nu; % Noise levels scaled to %
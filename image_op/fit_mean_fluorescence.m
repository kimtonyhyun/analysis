function M_norm = fit_mean_fluorescence(movie)
% Fit a sum of two exponentials to the mean flourescence
%
%   movie: movie matrix , [h x w x num_frames]
%
% 2015 01 31 Tony Hyun Kim (Revised: Hakan Inan, 15-Jan-4)
%

F = compute_mean_fluorescence(movie);
num_frames = length(F);

time = 1/20*(0:(num_frames-1));

f2 = fit(time',F','exp2');%fit sum of 2 exponentials
params = coeffvalues(f2);
F_fit = params(1)*exp(params(2)*time) + params(3)*exp(params(4)*time);

% Normalize movie
M_norm = zeros(size(movie), 'single');
for i = 1:num_frames
    M_norm(:,:,i) = single(movie(:,:,i))/F_fit(i);
end
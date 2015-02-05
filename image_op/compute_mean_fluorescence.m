function F = compute_mean_fluorescence(M)
% Compute the mean fluorescence of a movie M on a frame-by-frame basis.
%   M: [h x w x num_frames]
%
% 2015 01 31 Tony Hyun Kim

num_frames = size(M,3);
F = zeros(1,num_frames);
for i = 1:num_frames
    m = M(:,:,i);
    F(i) = mean(m(:));
end
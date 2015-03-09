function F = compute_fluorescence_stats(M)
% Compute fluorescence stats (min, mean, max) of a movie M on a 
% frame-by-frame basis.
%
% Inputs:
%   M: [height x width x num_frames] movie to be examined
%

num_frames = size(M,3);
F = zeros(num_frames, 3);
for i = 1:num_frames
    if (mod(i,2500)==0)
        fprintf('%s: Frames %d of %d examined...\n',...
            datestr(now), i, num_frames);
    end
    m = M(:,:,i);
    m = m(:);
    F(i,:) = [min(m) mean(m) max(m)];
end
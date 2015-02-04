function M_norm = norm_by_disk_filter(movie)
%Normalize every frame in the movie by a disk filtered version of itself
%
%   movie: movie matrix , [h x w x num_frames]
%
% 2015 01 31 Tony Hyun Kim (Revised: Hakan Inan, 15-Jan-4)
%

M_norm = zeros(size(movie), 'single');

% Apply spatial normalization
hDisk = fspecial('disk', 15);
m_f = imfilter(movie(:,:,1), hDisk, 'replicate');
m0 = mean(m_f(:));

for i = 1:num_frames
    if (mod(i,1000)==0)
        fprintf('  Frames %d of %d done\n', i, num_frames);
    end
    m_f = imfilter(movie(:,:,i), hDisk, 'replicate');
    m1 = mean(m_f(:));
    m_f = m0/m1*m_f;
    M_norm(:,:,i) = movie(:,:,i)./m_f;
end

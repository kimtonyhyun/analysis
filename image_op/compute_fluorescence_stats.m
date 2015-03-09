function F = compute_fluorescence_stats(M, varargin)
% Compute fluorescence stats (min, mean, max) of a movie M on a 
% frame-by-frame basis.
%
% Inputs:
%   M: [height x width x num_frames] movie to be examined
%

trim = 0;
if ~isempty(varargin)
    for k = 1:length(varargin)
        switch lower(varargin{k})
            case 'trim' % Trim borders
                trim = varargin{k+1};
        end
    end
end

num_frames = size(M,3);
F = zeros(num_frames, 3);
for frame_idx = 1:num_frames
    if (mod(frame_idx,2500)==0)
        fprintf('%s: Frames %d of %d examined...\n',...
            datestr(now), frame_idx, num_frames);
    end
    m = M((1+trim):(end-trim),...
          (1+trim):(end-trim),...
           frame_idx);
    m = m(:);
    F(frame_idx,:) = [min(m) mean(m) max(m)];
end
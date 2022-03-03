function Mb = load_behavior_movie_frames(vid, ind1, ind2)
% Loads a subset of the behavioral video frames.

num_frames = ind2 - ind1 + 1;

tic;
fprintf('%s: Loading %d frames from "%s"... ',...
    datestr(now), num_frames, vid.Name);
Mb = vid.read([ind1 ind2]);
Mb = squeeze(Mb(:,:,1,:));
t = toc;
fprintf('Done! (%.1f s; %.1f ms per frame)\n', t, t/num_frames*1e3);
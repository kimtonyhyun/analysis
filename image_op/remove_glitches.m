function M = remove_glitches(movie_in)
% Detect glitched frames by looking for anomalous minimum and maximum pixel
% values of each frame. Prompts the user for confirmation before replacing
% the movie frame.
%
% Usage:
%   M_gfix = remove_glitches('c11m3d19.hdf5');
%
% Todo:
%   - Interactive addition and removal of frames to be replaced

fprintf('%s: Loading movie...\n', datestr(now));
M = load_movie(movie_in);

fprintf('%s: Computing fluorescence stats...\n', datestr(now));
F = compute_fluorescence_stats(M);

z_thresh = 8;
fprintf('%s: Detect anomalous frames in the minimum trace...\n', datestr(now));
min_anomalous_frames = detect_anomalous_frames(F(:,1), z_thresh);

fprintf('%s: Detect anomalous frames in the max trace...\n', datestr(now));
max_anomalous_frames = detect_anomalous_frames(F(:,3), z_thresh);

num_frames = size(M,3);
plot(F(:,[1 3]), 'b.-');
xlim([1 num_frames]);
xlabel('Frame');
ylabel('Pixel values');
grid on;
hold on;
plot(min_anomalous_frames(:,1), min_anomalous_frames(:,2), 'ro');
plot(max_anomalous_frames(:,1), max_anomalous_frames(:,2), 'ro');
title('Red marker indicates frames to be replaced');

frames_to_replace = union(min_anomalous_frames(:,1),...
                          max_anomalous_frames(:,1));

fprintf('%s: Replace %d frames? (Press any key to continue) >> ',...
        datestr(now), length(frames_to_replace));
pause;

% Perform frame replacement
frames_to_replace = sort(frames_to_replace, 'descend')'; % Needs to be row vector
for frame = frames_to_replace
    if (frame == num_frames) % Last frame?
        source_frame = frame-1;
        while ismember(source_frame, frames_to_replace)
            source_frame = source_frame - 1;
        end
    else % Otherwise, just pull from the next frame
        source_frame = frame + 1;
    end
    M(:,:,frame) = M(:,:,source_frame);
    fprintf('  Frame %d replaced by frame %d\n',...
        frame, source_frame);
end

end % remove_glitches

function anomalous_frames = detect_anomalous_frames(trace, thresh_z)
% Detect anomalous frames of the one-dimensional signal `trace` by
% comparing the mean and standard deviation computed from the immediate
% neighborhood of the frame under test.

half_window = 30;
exclude_half_window = 5;

num_frames = length(trace);

anomalous_frames = [];
for i = 1:num_frames
    full_window = [i-half_window i+half_window];
    full_window(1) = max(1, full_window(1));
    full_window(2) = min(num_frames, full_window(2));
    
    exclude_window = (i-exclude_half_window):(i+exclude_half_window);
    nbhd = setdiff(full_window(1):full_window(2), exclude_window);
    nbhd_vals = trace(nbhd);
    
    mu = mean(nbhd_vals);
    sig = std(nbhd_vals);
    
    z = abs(trace(i) - mu)/sig;
    if (z > thresh_z)
        fprintf('  Frame %d has z-score of %.3f\n', i, z);
        anomalous_frames = [anomalous_frames; i trace(i)]; %#ok<AGROW>
    end
end

end % detect_anomalous_frames
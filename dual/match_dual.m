function match_dual(dir1, dir2)
% Examine frame count synchronicity between simultaneously acquired
% Miniscope recordings. Overwrite the TIF files in 'dir2' to match the
% frame counts determined from 'dir1'.
%
% Inputs:
%   - dir1: Directory containing TIF files from the "master" microscope.
%           The TIF files in this directory will NOT be altered.
%   - dir2: Directory containing TIF files from the "slave" microscope.
%           WARNING: TIF files in this directory may be edited so that the
%           frame counts align exactly to those in 'dir1'! Advice: Always
%           make a copy of the source data before applying this function!
%

frames1 = count_frames_in_tif(dir1);
[frames2, filenames2] = count_frames_in_tif(dir2);

num_files1 = length(frames1);
num_files2 = length(frames2);

assert(num_files1 == num_files2,...
       'Error, number of TIF files in "%s" (%d) does not match "%s" (%d)!',...
       dir1, num_files1, dir2, num_files2);
num_files = num_files1;
clear num_files1 num_files2;

frame_mismatch = frames2 - frames1;
max_mismatch = max(abs(frame_mismatch));

% Display results and prompt user to proceed
hist(frame_mismatch, -5:5);
grid on;
xlabel(sprintf('Frame mismatch ("%s"-"%s")', dir2, dir1));
ylabel('Number of instances');
title(sprintf('Num file pairs=%d; Max mismatch=%d frames', num_files, max_mismatch));

fprintf('match_tifs_dual: WARNING! This function will modify TIF files in place!\n');
input('match_tifs_dual: Press enter to proceed >> ');

% Begin editing of TIF files
for k = 1:num_files
    filename2 = filenames2{k};
    mismatch = frame_mismatch(k);
    if (mismatch ~= 0)
        fprintf('%s: Editing file "%s" with %d extra frames...\n',...
            datestr(now), filename2, mismatch);
        M_orig = load_movie(filename2);
        
        if (mismatch > 0)
            % File2 has extra frames. In this case, we chop off the extra
            % frames from the beginning of the trial.
            M_repl = M_orig(:,:,(1+mismatch):end);
        elseif (mismatch < 0)
            % File2 has fewer frames. In this case, we duplicate extra
            % frames at the end of the trial.
            frame = M_orig(:,:,end);
            M_extra = repmat(frame, 1, 1, abs(mismatch));
            M_repl = cat(3, M_orig, M_extra);
        end
        
        % Replace TIF file!
        delete(filename2);
        save_movie_to_tif(M_repl, filename2);
    end
end
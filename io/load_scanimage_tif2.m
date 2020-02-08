function M = load_scanimage_tif2(filename, varargin)
% Used for loading ScanImage TIF files on Linux. Tested on Matlab R2019a on
% Ubuntu 18.04
%

info = imfinfo(filename);
num_total_frames = length(info);
fprintf('  File "%s" has %d frames\n', filename, num_total_frames);

frames_to_load = 1:num_total_frames;
for k = 1:length(varargin)
    if ischar(varargin{k})
        switch lower(varargin{k})
            case 'odd'
                fprintf('  Loading ODD frames only!\n');
                frames_to_load = 1:2:num_total_frames;
            case 'even'
                fprintf('  Loading EVEN frames only!\n');
                frames_to_load = 2:2:num_total_frames;
            case 'frames'
                frames_to_load = varargin{k+1};
        end
    end
end
num_frames = length(frames_to_load);

width = info(1).Width;
height = info(1).Height;
M = zeros(height, width, num_frames, 'int16');

q = 1;
for k = frames_to_load
    if (mod(q,1000)==0)
        fprintf('  %s: Frames %d / %d loaded\n', datestr(now), q, num_frames);
    end
    M(:,:,q) = imread(filename, k, 'info', info);
    q = q + 1;
end
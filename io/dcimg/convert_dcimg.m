function movie_out = convert_dcimg(movie_in, movie_out)
% Convert Hamamatsu DCIMG into a HDF5 file. TODO: Is converting one frame
% at a time truly the only way to access DCIMG data?
%

if ~exist('movie_out', 'var')
    [~, movie_out] = fileparts(movie_in);
    movie_out = sprintf('%s.hdf5', movie_out);
end
fprintf('convert_dcimg: Output movie will be saved as "%s"\n', movie_out);

% Analyze DCIMG file
%------------------------------------------------------------
data_type = 'uint16';
[A, num_frames] = dcimgmatlab(int32(0), source);
[height, width] = size(A);

height = double(height);
width = double(width);
num_frames = double(num_frames);

% Prep output
%------------------------------------------------------------
dataset_name = '/Data/Images';

h5create(movie_out, dataset_name, [height width num_frames],...
    'DataType', data_type,...
    'ChunkSize', [height width 1]);

for k = 1:num_frames
    if (mod(k,1000)==0)
        fprintf('  %s: At frames %d / %d\n', datestr(now), k, num_frames);
        
        % DCIMG MEX accumulates in memory. Need to clear it out
        clear mex; %#ok<*CLMEX>
    end
    
    A = dcimgmatlab(int32(k-1), source); % 0-based indexing
    h5write(movie_out, dataset_name, A, [1 1 k], [height width 1]);
end

% Create standard "/Params" directory
h5create(movie_out, '/Params/NumFrames', 1);
h5write(movie_out, '/Params/NumFrames', num_frames);

clear mex; % Final cleanup
fprintf('%s: Done!\n', datestr(now));
function movie = load_tif_to_memory(source)

% Usage: M = load_tif_to_memory('mouse1_day4_sp2_mc_cr.tif');

info = imfinfo(source);
num_frames = length(info);
width  = info(1).Width;
height = info(1).Height;

movie = zeros(height, width, num_frames, 'single');

% Load into memory
t = Tiff(source, 'r');
for k = 1:num_frames
    if (mod(k,1000)==0)
        fprintf('  Frames %d / %d loaded\n', k, num_frames);
    end
    t.setDirectory(k);
    movie(:,:,k) = t.read();
end
t.close();
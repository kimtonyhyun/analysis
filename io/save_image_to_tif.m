function save_image_to_tif(A, filename)
% Export TIFF images in single format, for loading into Photoshop
%
% Note: After loading image into Photoshop, to see the image content:
%   Image -> Adjustments -> HDR Toning
%   Try 'Equalize histogram' or 'Local adaptation' methods
%
% TODO:
%   - Other types of outputs, e.g. color-fusion
%
t = Tiff(filename, 'w');
[h, w] = size(A);
t.setTag('ImageWidth', w);
t.setTag('ImageLength', h);
t.setTag('Photometric', Tiff.Photometric.MinIsBlack);
t.setTag('Compression', Tiff.Compression.None);
t.setTag('BitsPerSample', 32);
t.setTag('SamplesPerPixel', 1);
t.setTag('PlanarConfiguration', Tiff.PlanarConfiguration.Chunky); 
t.setTag('SampleFormat', Tiff.SampleFormat.IEEEFP);
write(t,A);

close(t);
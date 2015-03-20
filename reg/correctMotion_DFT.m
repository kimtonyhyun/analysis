function correctMotion_DFT(hdf5_dir,hdf5_in_name,hdf5_out_name)

dataset = '/Data/Images';

h5att = h5info(fullfile(hdf5_dir,hdf5_in_name),dataset);
imageStackSize = h5att.Dataspace.Size;
rows = imageStackSize(1);
cols = imageStackSize(2);
frames = imageStackSize(3);
baseImage = h5read(fullfile(hdf5_dir,hdf5_in_name),dataset,[1,1,1],[rows,cols,1]);

ssm_radius = 20;
asm_radius = 5;

hDisk  = fspecial('disk', ssm_radius);
hDisk2 = fspecial('disk', asm_radius);
transform = @(A) mosaic_transform(A, hDisk, hDisk2);

im_ref = transform(baseImage);

% Specify ROI
%------------------------------------------------------------
h_roi = figure;
imagesc(im_ref); axis image; colormap gray;
title('Select ROI');
h_rect = imrect;
mask_rect = round(getPosition(h_rect));
cropBase = im_ref(mask_rect(2):mask_rect(2)+mask_rect(4),mask_rect(1):mask_rect(1)+mask_rect(3));
close(h_roi);

baseFFT = fft2(cropBase);

chunkSize = [rows,cols,1];
h5create(fullfile(hdf5_dir,hdf5_out_name),dataset,[Inf,Inf,Inf],'ChunkSize',chunkSize,'Datatype','single');
h5write(fullfile(hdf5_dir,hdf5_out_name),dataset,baseImage,[1,1,1],[rows,cols,1]);

h = waitbar(0,'Motion Correction in Progress...');
for f=2:frames
    waitbar(f/frames);
    image = h5read(fullfile(hdf5_dir,hdf5_in_name),dataset,[1,1,f],[rows,cols,1]);
    cropImage = image(mask_rect(2):mask_rect(2)+mask_rect(4),mask_rect(1):mask_rect(1)+mask_rect(3));
    tImage = transform(cropImage);
    imageFFT = fft2(tImage);
    [output, result] = dftregistration(baseFFT,imageFFT,100); 
    
    %%% Apply shift parameters to the non-blurred/cropped image
    iFFT = fft2(image);
    [nr,nc]=size(iFFT);
    Nr = ifftshift([-fix(nr/2):ceil(nr/2)-1]);
    Nc = ifftshift([-fix(nc/2):ceil(nc/2)-1]);
    [Nc,Nr] = meshgrid(Nc,Nr);
    iFFT = iFFT.*exp(i*2*pi*(-output(3)*Nr/nr-output(4)*Nc/nc));
    regImageFFT = iFFT*exp(i*output(2));
    regImage = abs(ifft2(regImageFFT));
        
    h5write(fullfile(hdf5_dir,hdf5_out_name),dataset,regImage,[1,1,f],[rows,cols,1]);
end
close(h);
end

function A_tr = mosaic_transform(A, ssm_filter, asm_filter)
    A_tr = A - imfilter(A, ssm_filter, 'replicate');
    A_tr = imfilter(A_tr, asm_filter);
end

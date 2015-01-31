function M_reg = run_turboreg(M)

[h, w, num_frames] = size(M);

% Select reference frame
%------------------------------------------------------------
fprintf('%s: Specify parameters for TurboReg\n', datestr(now));
ref_frame_idx = get_integer_input('  Select reference frame >> ',...
                                  1, [1 num_frames]);
fprintf('    -> Frame %d selected as registration reference\n', ref_frame_idx);
im_ref = M(:,:,ref_frame_idx);

% Select "Mosaic filter" parameters
%------------------------------------------------------------
ssm_radius = get_integer_input('  Specify radius for "subtract spatial mean" >> ', 20, [5 50]);
fprintf('    -> Subtract spatial mean radius set to %d\n', ssm_radius);
asm_radius = get_integer_input('  Specify radius for "apply spatial mean" >> ', 5, [1 20]);
fprintf('    -> Apply spatial mean radius set to %d\n', asm_radius);

hDisk  = fspecial('disk', ssm_radius);
hDisk2 = fspecial('disk', asm_radius);
transform = @(A) mosaic_transform(A, hDisk, hDisk2);

im_ref = transform(im_ref);

% Specify ROI
%------------------------------------------------------------
h_roi = figure;
imagesc(im_ref); axis image; colormap gray;
title('Select ROI');
h_poly = impoly;
mask_poly = getPosition(h_poly);
mask_ref = single(poly2mask(mask_poly(:,1), mask_poly(:,2), h, w));
close(h_roi);

% TurboReg parameters
%------------------------------------------------------------
options.rotation_enable = true;
options.mingain = 0.0; % Max accuracy
options.levels = calculate_pyramid_depth(min(w, h));

fprintf('%s: Began TurboReg registration with %d pyramid levels...\n',...
    datestr(now), options.levels);
M_reg = zeros(size(M), class(M));
tic;
parfor i = 1:num_frames
    im_reg = transform(M(:,:,i));
    im_coreg = M(:,:,i);
    M_reg(:,:,i) = turbocoreg(im_ref, mask_ref, im_reg, im_coreg, options);
end
t = toc;
fprintf('%s: Finished in %.1f sec (%.1f sec per frame)\n',...
    datestr(now), t, t/num_frames);

end

function val = get_integer_input(prompt, default, limits)
    while (1)
        val = str2double(input(prompt, 's'));
        if (~isnan(val))
            if ((limits(1) <= val) && (val <= limits(2)))
                val = floor(val);
                break;
            end
        else
            val = default;
            break;
        end
    end
end

function A_tr = mosaic_transform(A, ssm_filter, asm_filter)
    A_tr = A - imfilter(A, ssm_filter, 'replicate');
    A_tr = imfilter(A_tr, asm_filter);
end

function depth = calculate_pyramid_depth(len)
    min_size = 45;
    depth = 0;
    while (min_size <= len)
        len = len/2;
        depth = depth + 1;
    end
end
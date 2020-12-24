function [colorMasks, chartMask] = getColorCheckerMasks(imageXYZ, allowRoiSelect, configuration)

% Use the provided estimate of ROI center coordinates.
switch configuration 
    case 'singleLightFar'
        xx = [ 33, 70, 103, 140, 175, 212];
        yy = [ 39, 76, 109, 147, 182];
        r = 15;
        isRotated = false; 
    case 'singleLightClose'
        xx = [47, 89, 139, 182, 233, 279];
        yy = [38, 93, 140, 184, 229];
        r = 20;
        isRotated = true; 
    case 'doubleLightClose'
        xx = [47, 89, 139, 182, 233, 279];
        yy = [38, 93, 140, 184, 229];
        r = 20;
        isRotated = true; 
    otherwise
        xx = [47, 89, 139, 182, 233, 279];
        yy = [38, 93, 140, 184, 229];
        isRotated = true;
        wargning('Roi center coordinates need re-evaluation');
end 

saveFilename = mkNewDir(getSetting('matdir'), configuration, strcat(strrep(configuration, '.h5', ''), '_others.mat'));

if nargin < 2
    allowRoiSelect = false;
end

imageY = imageXYZ(:, :, 2);

warning('off', 'images:initSize:adjustingMag')
figure(1);
imagesc(imageY);
title('Tristimulus Y image');

if allowRoiSelect
    % Click to add vertices, then right-click and select "Create Mask" to return.
    title('Draw a polygon around the chart')
    chartMask = roipoly; %When you are finished positioning and sizing the polygon, create the mask by double-clicking, or by right-clicking inside the region and selecting Create mask from the context menu.
    chartMask = imdilate(chartMask, ones(7));
    save(saveFilename, 'chartMask');
%     load(strcat(saveFilename), 'chartMask');

    croppedXYZ = imageXYZ(any(chartMask, 2), any(chartMask, 1), :);
    croppedY = croppedXYZ(:, :, 2);

    A = croppedXYZ;
    figure(1);
    imagesc(croppedY);
    title('Cropped ROI for color chart');

    [x, y] = meshgrid(xx, yy);

    if isRotated
        tmp = yy;
        yy = xx;
        xx = tmp;
        [x, y] = ndgrid(xx, yy);
    end

    x = x';
    y = y';
    x = x(:);
    y = y(:);

    x = round(x);
    y = round(y);

    mask = false(size(A, 1), size(A, 2));
    colorMasks = false(size(A, 1), size(A, 2), length(x));
    for k = 1:length(x)
        mask(y(k)-r:y(k)+r, x(k)-r:x(k)+r) = true;

        colorMask = false(size(A, 1), size(A, 2));
        colorMask(y(k)-r:y(k)+r, x(k)-r:x(k)+r) = true;
        colorMasks(:, :, k) = colorMask;
        imshow(colorMask);
    end

    mask_eroded = imerode(mask, strel('disk', 5));

    %mask_clipped = (A == intmax(class(A))) | (A == intmin(class(A)));
    mask_clipped = (A == intmax('uint8')) | (A == intmin('uint8'));
    mask_clipped = mask_clipped(:, :, 1) | mask_clipped(:, :, 2) | mask_clipped(:, :, 3);

    mask_patches = mask_eroded & ~mask_clipped;

    baseImage = ind2rgb(uint8(croppedY), parula(255));
    A_patches = imoverlay(baseImage, mask_patches);

    figure(2);
    imshow(A_patches);
    title('The selected pixels are highlighted in yellow');
    save(saveFilename, 'colorMasks', 'chartMask');
else 
    load(saveFilename, 'colorMasks', 'chartMask');
end 


end
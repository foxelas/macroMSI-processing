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
    case 'capture_average_comparison'
    case 'fusion_comparison'
        xx = [41, 93, 147, 195, 248, 295];
        yy = [44, 94, 144, 195, 245];
        isRotated = true;
        r = 1; % 1, 20
    case 'polarizing_effect_on_tissue'
        xx = [41, 93, 147, 195, 248, 295];
        yy = [44, 94, 144, 195, 245];
        isRotated = true;
        r = 1;
    otherwise
        xx = [41, 93, 147, 195, 248, 295];
        yy = [44, 94, 144, 195, 245];
        isRotated = true;
        r = 20;
        warning('Roi center coordinates need re-evaluation');
end 

saveFilename = mkNewDir(getSetting('matdir'), getSetting('experiment'), strcat(strrep(configuration, '.h5', ''), '_others.mat'));

if nargin < 2
    allowRoiSelect = false;
end

has3Channels = ndims(imageXYZ) == 3;
if has3Channels
    imageY = imageXYZ(:, :, 2);
else 
    imageY = imageXYZ;
end

warning('off', 'images:initSize:adjustingMag')
figure(1);
imagesc(imageY);
title('Tristimulus Y image');

if allowRoiSelect
    if ~exist(saveFilename, 'file')
        % Click to add vertices, then right-click and select "Create Mask" to return.
        title('Draw a polygon around the chart')
        chartMask = roipoly; %When you are finished positioning and sizing the polygon, create the mask by double-clicking, or by right-clicking inside the region and selecting Create mask from the context menu.
        chartMask = imdilate(chartMask, ones(7));
        save(saveFilename, 'chartMask');
        fprintf('Saved file for color chart mask.\n');
    else 
        load(saveFilename, 'chartMask');
        fprintf('Loaded file for color chart mask.\n');
    end 
%     load(strcat(saveFilename), 'chartMask');

    croppedY = imageY(any(chartMask, 2), any(chartMask, 1));
    if has3Channels
        croppedXYZ =  imageY(any(chartMask, 2), any(chartMask, 1), :);
        A = croppedXYZ;
    else
        A = croppedY;
    end 
%     figure(1);
%     imagesc(croppedY);
%     title('Cropped ROI for color chart');

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

    figure(1);
    mask = false(size(A, 1), size(A, 2));
    colorMasks = false(size(A, 1), size(A, 2), length(x));
    for k = 1:length(x)
        mask(y(k)-r:y(k)+r, x(k)-r:x(k)+r) = true;

        colorMask = false(size(A, 1), size(A, 2));
        colorMask(y(k)-r:y(k)+r, x(k)-r:x(k)+r) = true;
        colorMasks(:, :, k) = colorMask;
        imshow(colorMask);
    end

    if r == 1 
        mask_eroded = mask; 
    else 
        mask_eroded = imerode(mask, strel('disk', 5));
    end 
    
    %mask_clipped = (A == intmax(class(A))) | (A == intmin(class(A)));
    mask_clipped = (A == intmax('uint8')) | (A == intmin('uint8'));
    if has3Channels
        mask_clipped = mask_clipped(:, :, 1) | mask_clipped(:, :, 2) | mask_clipped(:, :, 3);
    end
    mask_patches = mask_eroded & ~mask_clipped;

    if has3Channels
        baseImage = ind2rgb(uint8(croppedY), parula(255));
        A_patches = imoverlay(baseImage, mask_patches);
    else
        A_patches = imoverlay(rescale(croppedY), mask_patches);
    end

    figure(5);
    imagesc(A_patches);
    title('The selected pixels are highlighted in yellow');
    save(saveFilename, 'colorMasks', 'chartMask');
else 
    load(saveFilename, 'colorMasks', 'chartMask');
end 


end
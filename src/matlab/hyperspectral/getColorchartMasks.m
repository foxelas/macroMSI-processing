function [colorMasks, chartMask] = getColorchartMasks(baseImage, allowRoiSelection, configuration)
%GETCOLOCHECKERMASKS returns masks for the colorchart
%
%   [colorMasks, chartMask] = getColorchartMasks(baseImage,
%     allowRoiSelection, configuration) returns colorchart masks for the
%     specific baseImage and configuration
%

if nargin < 2
    allowRoiSelection = false;
end

if nargin < 3
    configuration = 'unknown';
end

%% Settings for locations of colorchart patch centers
% Use the provided estimate of ROI center coordinates.
switch configuration
    case 'singleLightFar'
        xx = [33, 70, 103, 140, 175, 212];
        yy = [39, 76, 109, 147, 182];
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

if isempty(getSetting('targetPosition'))
    suffix = '';
else
    suffix = strcat('_', getSetting('targetPosition'));
end

saveFilename = mkNewDir(getSetting('matdir'), getSetting('experiment'), ...
    strcat(strrep(configuration, '.h5', ''), '_others', suffix, '.mat'));

has3Channels = ndims(baseImage) == 3;
if has3Channels
    baseImage = baseImage(:, :, 2);
end

warning('off', 'images:initSize:adjustingMag')
figure(1);
imagesc(baseImage);
title('Base Image');

if allowRoiSelection
    if ~exist(saveFilename, 'file')
        % Click to add vertices, then right-click and select "Create Mask" to return.
        title('Draw a polygon around the chart')
        chartMask = roipoly;
        %When you are finished positioning and sizing the polygon, create the mask by double-clicking, or by right-clicking inside the region and selecting Create mask from the context menu.
        chartMask = imdilate(chartMask, ones(7));
        save(saveFilename, 'chartMask');
        fprintf('Saved file for color chart mask.\n');
    else
        load(saveFilename, 'chartMask');
        fprintf('Loaded file for color chart mask.\n');
    end
    
    if has3Channels
        croppedImage = baseImage(any(chartMask, 2), any(chartMask, 1), :);
    else
        croppedImage = baseImage(any(chartMask, 2), any(chartMask, 1));
    end
    A = croppedImage;
    
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
    
    figure(2);
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
        baseImage = ind2rgb(uint8(croppedImage), parula(255));
        A_patches = imoverlay(baseImage, mask_patches);
    else
        A_patches = imoverlay(rescale(croppedImage), mask_patches);
    end
    
    figure(3);
    imagesc(A_patches);
    title('The selected pixels are highlighted in yellow');
    save(saveFilename, 'colorMasks', 'chartMask');
else
    load(saveFilename, 'colorMasks', 'chartMask');
end


end
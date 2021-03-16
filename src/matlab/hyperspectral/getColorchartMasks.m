function [colorMasks, chartMask] = getColorchartMasks(baseImage, allowRoiSelection, filename)
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
    filename = 'unknown';
    warning('Unknown filename, saving in unknown.mat');
end

%% Settings for locations of colorchart patch centers
% Use the provided estimate of ROI center coordinates.
configuration = getSetting('experiment');
isRotated = getSetting('isRotated');
switch configuration
    case 'singleLightFar'
        xx = [40, 70, 103, 140, 175, 212];
        yy = [39, 76, 109, 147, 182];
        r = 15;
        isRotated = false;
    case 'testStomach'
        xx = [28, 53, 77, 105, 129, 153]-2;
        yy = [25, 52, 76, 101, 134]-2;
        r = 0;
%     case 'singleLightClose'
%         xx = [47, 89, 139, 182, 233, 279];
%         yy = [38, 93, 140, 184, 229];
%         r = 20;
%         isRotated = true;
%     case 'doubleLightClose'
%         xx = [47, 89, 139, 182, 233, 279];
%         yy = [38, 93, 140, 184, 229];
%         r = 20;
%         isRotated = true;
%     case 'capture_average_comparison'
%     case 'fusion_comparison'
%         xx = [41, 93, 147, 195, 248, 295];
%         yy = [44, 94, 144, 195, 245];
%         isRotated = true;
%         r = 1; % 1, 20
%     case 'polarizing_effect_on_tissue'
%         xx = [41, 93, 147, 195, 248, 295];
%         yy = [44, 94, 144, 195, 245];
%         isRotated = true;
%         r = 1;
    otherwise
%         xx = [41, 93, 147, 195, 248, 295];
%         yy = [44, 94, 144, 195, 245];
        xx = [47, 89, 139, 182, 233, 279];
        yy = [44, 93, 140, 184, 229];
        r = 0;
end

saveFilename = mkNewDir(getSetting('matdir'), 'cororchartMasks', ...
    strcat('colorchartMask_', filename,'.mat'));

warning('off', 'images:initSize:adjustingMag')
fig1 = figure(1);
imshow(baseImage);
title('Base Image');

if allowRoiSelection
    setSetting('plotName', mkNewDir(getSetting('savedir'), getSetting('saveFolder'), filename));
    savePlot(fig1);
    
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
    
 
    A = baseImage(any(chartMask, 2), any(chartMask, 1), :);
    
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
    
    clf(fig1);
    mask = false(size(A, 1), size(A, 2));
    colorMasks = false(size(A, 1), size(A, 2), length(x));
    for k = 1:length(x)
        mask(y(k)-r:y(k)+r, x(k)-r:x(k)+r) = true;
        
        colorMask = false(size(A, 1), size(A, 2));
        colorMask(y(k)-r:y(k)+r, x(k)-r:x(k)+r) = true;
        colorMasks(:, :, k) = colorMask;
        imshow(colorMask);
    end
    
    if r < 2
        mask_eroded = mask;
    else
        mask_eroded = imerode(mask, strel('disk', 5));
    end
    
%     %mask_clipped = (A == intmax(class(A))) | (A == intmin(class(A)));
%     mask_clipped = (A == intmax('uint8')) | (A == intmin('uint8'));
%     mask_clipped = mask_clipped(:, :, 1) | mask_clipped(:, :, 2) | mask_clipped(:, :, 3);
%     mask_patches = mask_eroded & ~mask_clipped;
    mask_patches = mask_eroded;
    A_patches = imoverlay(A, mask_patches);
    
    clf(fig1);
    imshow(A_patches);
    title('The selected pixels are highlighted in yellow');
    setSetting('plotName', mkNewDir(getSetting('savedir'), getSetting('saveFolder'), strcat(filename, '_patchMask')));
    savePlot(fig1);
    save(saveFilename, 'colorMasks', 'chartMask');
else
    load(saveFilename, 'colorMasks', 'chartMask');
end


end
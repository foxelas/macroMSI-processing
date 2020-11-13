function [colorMasks, chartMask] = getColorCheckerMasks(imageXYZ, fileName, allowRoiSelect)

isRotated = true;
if nargin < 3
    allowRoiSelect = false;
end

imageY = imageXYZ(:, :, 2);

warning('off', 'images:initSize:adjustingMag')
figure(1);
imagesc(imageY);
title('Tristimulus Y image');

if allowRoiSelect
    % Use roipoly to create a mask from a polygon drawn manually.
    % Click to add vertices, then right-click and select "Create Mask" to return.
    title('Draw a polygon around the chart')
    chartMask = roipoly; %When you are finished positioning and sizing the polygon, create the mask by double-clicking, or by right-clicking inside the region and selecting Create mask from the context menu.
    chartMask = imdilate(chartMask, ones(7));
    save(strcat(fileName, '.mat'), 'chartMask');
end

load(strcat(fileName, '.mat'), 'chartMask');

croppedXYZ = imageXYZ(any(chartMask, 2), any(chartMask, 1), :);
croppedY = croppedXYZ(:, :, 2);

A = croppedXYZ;
figure(1);
imagesc(croppedY);
title('Cropped ROI for color chart');

% Use the provided estimate of ROI center coordinates.
xx = [36, 75, 106, 143, 178, 215] - 3;
yy = [45, 82, 115, 153, 188] - 7;
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
%r = mean(diff(x)) / 2 * 0.60;
%r = floor(r);
r = 10;

mask = false(size(A, 1), size(A, 2));
colorMasks = false(size(A, 1), size(A, 2), length(x));
for k = 1:length(x)
    mask(y(k)-r:y(k)+r, x(k)-r:x(k)+r) = true;

    colorMask = false(size(A, 1), size(A, 2));
    colorMask(y(k)-r:y(k)+r, x(k)-r:x(k)+r) = true;
    colorMasks(:, :, k) = colorMask;
end

mask_eroded = imerode(mask, strel('disk', 5));
completeMask = mask_eroded;

%mask_clipped = (A == intmax(class(A))) | (A == intmin(class(A)));
mask_clipped = (A == intmax('uint8')) | (A == intmin('uint8'));
mask_clipped = mask_clipped(:, :, 1) | mask_clipped(:, :, 2) | mask_clipped(:, :, 3);

mask_patches = mask_eroded & ~mask_clipped;

baseImage = ind2rgb(uint8(croppedY), parula(255));
A_patches = imoverlay(baseImage, mask_patches);

figure(2);
imshow(A_patches);
title('The selected pixels are highlighted in yellow')
end
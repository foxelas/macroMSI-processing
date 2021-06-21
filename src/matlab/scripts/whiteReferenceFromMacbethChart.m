close all;
locationDir = 'D:\temp\mspi\saitamav2\white\unfixed\right';
whiteFilename = 'P2270082.tif';
msiFilenames = {'P2270083.tif', 'P2270084.tif', 'P2270085.tif', 'P2270086.tif', 'P2270087.tif', 'P2270088.tif', 'P2270089.tif'};
allFilenames = [whiteFilename, msiFilenames];
allFilenames = arrayfun(@(x) fullfile(locationDir, x), allFilenames);

%% Test various color correction techniques on white light image
A = imread(fullfile(locationDir, whiteFilename));
warning('off', 'images:initSize:adjustingMag')
load('parameters\mask_chart.mat', 'mask_chart');

% Use the provided estimate of ROI center coordinates.
x = 1262;
y = 2498;

x = round(x);
y = round(y);
r = 145;

mask = false(size(A, 1), size(A, 2));
mask(y-r:y+r, x-r:x+r) = true;

mask_eroded = imerode(mask, strel('disk', 5));

mask_clipped = (A == intmax(class(A))) | (A == intmin(class(A)));
mask_clipped = mask_clipped(:, :, 1) | mask_clipped(:, :, 2) | mask_clipped(:, :, 3);

mask_patches = mask_eroded & ~mask_clipped;

A_patches = imoverlay(A, mask_patches);

figure(1); imshow(A_patches)
title('The selected pixels are highlighted in yellow');

height = length(y-r:y+r);
width = length(x-r:x+r);
coordinates = [x - r, y - r,];

[referencePatchMSI, ~, ~, segmentMask, segmentMaskI] = readMSI(allFilenames, coordinates, width, height);
columnMSI = reshape(referencePatchMSI, size(referencePatchMSI, 1), size(referencePatchMSI, 2)*size(referencePatchMSI, 3), size(referencePatchMSI, 4));
referenceWhite = mean(columnMSI, 2);
save(getSetting('whiteReferenceMacbeth'), 'referenceWhite', 'referencePatchMSI');
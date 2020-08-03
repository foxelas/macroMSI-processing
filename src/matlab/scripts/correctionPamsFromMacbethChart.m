% close all;
% locationDir = 'D:\temp\mspi\saitamav2\white\unfixed\right';
% setSetting('saveImages', false); 
% whiteFilename = 'P2270082.tif';
% msiFilenames = {'P2270083.tif', 'P2270084.tif', 'P2270085.tif', 'P2270086.tif', 'P2270087.tif', 'P2270088.tif', 'P2270089.tif'};
% allFilenames = [whiteFilename, msiFilenames];
% allFilenames = arrayfun(@(x) fullfile(locationDir, x), allFilenames);
% 
% %% Test various color correction techniques on white light image
% A = imread(fullfile(locationDir, whiteFilename));
% warning('off', 'images:initSize:adjustingMag')
% load('parameters\mask_chart.mat', 'mask_chart');
% 
% % Use the provided estimate of ROI center coordinates.
% r = 145;
% 
% %% for red 
% x = 2262;
% y = 1992;
% [x, y] = showSelectedPatch(A, x, y, r);
% [referencePatchMSI, referenceAvg] = readImagePatch(allFilenames, x, y, r);
% redCoeff = im2double(im2uint8(referenceAvg(7, :, 1))) / im2double(uint8(175));
% 
% %% for green 
% x = 1765;
% y = 1986;
% [x, y] = showSelectedPatch(A, x, y, r);
% [referencePatchMSI, referenceAvg] = readImagePatch(allFilenames, x, y, r);
% greenCoeff = im2double(im2uint8(referenceAvg(5, :, 2))) / im2double(uint8(148));


%% for blue 
x = 1272;
y = 1974;
[x, y] = showSelectedPatch(A, x, y, r);
[referencePatchMSI, referenceAvg] = readImagePatch(allFilenames, x, y, r);
blueCoeff = im2double(im2uint8(referenceAvg(1, :, 3))) / im2double(uint8(150));


save(getSetting('channelCorrection'), 'redCoeff', 'greenCoeff', 'blueCoeff');

function [x, y] = showSelectedPatch(A, x, y, r)
    x = round(x);
    y = round(y);

    mask = false(size(A, 1), size(A, 2));
    mask(y-r:y+r, x-r:x+r) = true;

    mask_eroded = imerode(mask, strel('disk', 5));

    mask_clipped = (A == intmax(class(A))) | (A == intmin(class(A)));
    mask_clipped = mask_clipped(:, :, 1) | mask_clipped(:, :, 2) | mask_clipped(:, :, 3);

    mask_patches = mask_eroded & ~mask_clipped;

    A_patches = imoverlay(A, mask_patches);

    figure(1); imshow(A_patches)
    title('The selected pixels are highlighted in yellow');
end 

function  [referencePatchMSI, referenceAvg] = readImagePatch(allFilenames, x, y, r) 
    height = length(y-r:y+r); 
    width = length(x-r:x+r);
    coordinates = [x-r, y-r,];
    [referencePatchMSI, ~, ~, ~, ~] = readMSI(allFilenames, coordinates, width, height); 
    columnMSI = reshape(referencePatchMSI, size(referencePatchMSI, 1), size(referencePatchMSI, 2) * size(referencePatchMSI, 3), size(referencePatchMSI, 4) );
    referenceAvg = mean(columnMSI, 2); 
end 
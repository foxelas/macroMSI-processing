disablePrepare = true;
if ~disablePrepare
    prepareData;
end

disable = true; 
allowRoiSelection = false;

filename = getFilename(configuration, 'colorchart');
patchesStart = 1;
patchesEnd = length(spectraColorOrder);
selectedPatches = [patchesStart:patchesEnd];

if ~disable 
    [spectralData, imageXYZ, wavelengths] = loadH5Data(filename, configuration);
    [colorMasks, chartMask] = getColorCheckerMasks(imageXYZ, allowRoiSelection, configuration);
    croppedSpectral = getHSIdata(spectralData, chartMask);
    [actualSpectralVals] = getPatchValues(croppedSpectral, colorMasks);
end         

[reorderedSpectralVals, lineNames] = ReorderSpectra(actualSpectralVals, chartColorOrder, spectraColorOrder, wavelengths, expectedWavelengths);
[reorderedSpectralValsRaw, ~] = ReorderSpectra(actualSpectralVals, chartColorOrder, spectraColorOrder, wavelengths, expectedWavelengths);

%% Standard (expected) color patch spectra
plotFunWrapper(4, @plotColorChartSpectra, expectedWavelengths, expectedSpectra(selectedPatches, :), lineNames, {saveFolder, 'expected'});
plotFunWrapper(6, @plotColorChartSpectra, expectedWavelengths, reorderedSpectralVals(selectedPatches, :), lineNames, {saveFolder, 'measured'});

gofs = applyFuncOnRows(reorderedSpectralVals(selectedPatches, :), expectedSpectra(selectedPatches, :), @goodnessOfFit)
nmses = applyFuncOnRows(reorderedSpectralVals(selectedPatches, :), expectedSpectra(selectedPatches, :), @nmse)

save(fullfile(getSetting('savedir'), getSetting('saveFolder'), 'last_run.mat'), '-v7.3');

    %     %% Get location of color chart
    %     croppedXYZ = imageXYZ(any(chartMask,2), any(chartMask, 1), :);
    %     croppedXYZnorm = croppedXYZ ./ sum(sum(croppedXYZ(:,:,2)));
    %     rgb = xyz2rgb(croppedXYZnorm,'ColorSpace','srgb');
    %
    %     croppedLab = xyz2lab(imageXYZ(any(chartMask,2), any(chartMask, 1), :));
    %     %% Get average RGB of each color
    %
    %     %% Compare measured RGB with real RGB
    %     % for RGB comparison
    %     %euclidDiff = sqrt((L1-L2).^2 + (a1-a2).^2 + (b1-b2).^2);
    %     dE = imcolordiff(I1,I2, 'Standard', 'CIEDE2000');
    %
    %     %Display the color difference as an image.
    %     %Scale the display range to match the range of pixel values in dE.
    %     %Bright regions indicate the greatest color difference and correspond with the pink regions of tissue.
    %
    %     imshow(dE,[]);
    %
    %
    %     %% Get average Spectrum of each color
    %
    %     %% Compare measured Spectrum with real Spectrum


function [reorderedSpectra, labels] = ReorderSpectra(target, chartColorOrder, spectraColorOrder, wavelengths, spectralWavelengths)

[~, idx] = ismember(spectraColorOrder, chartColorOrder);
idx = nonzeros(idx); 
[~, idx2] = ismember(spectralWavelengths', wavelengths);

targetDecim = target(:, idx2);
reorderedSpectra = targetDecim(idx, :);

labels = spectraColorOrder;
end

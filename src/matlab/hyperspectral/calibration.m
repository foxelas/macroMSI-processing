toRun = false;
if toRun
    prepareData;
end

iStart = getIndexValue('macbethFile', dataDate);
patchesStart = 1;
patchesEnd = 24;
selectedPatches = [patchesStart:patchesEnd];

whites = zeros(length(wavelengths), length(datafiles));
for i = iStart:length(datafiles)
    filename = datafiles(i).name;
    
    if toRun 
        [spectralData, imageXYZ, wavelengths] = loadH5Data(filename);
        [colorMasks, chartMask] = getColorCheckerMasks(imageXYZ, filename, true);
        croppedSpectral = spectralData(any(chartMask, 2), any(chartMask, 1), :);
        [actualSpectralVals] = getPatchValues(croppedSpectral, colorMasks);
    end 

    [reorderedSpectralVals, lineNames] = ReorderSpectra(actualSpectralVals, chartColorOrder, spectraColorOrder, wavelengths, expectedWavelengths);
    [reorderedSpectralValsRaw, ~] = ReorderSpectra(actualSpectralVals, chartColorOrder, spectraColorOrder, wavelengths, expectedWavelengths, true);

    %% Standard (expected) color patch spectra
    plotFunWrapper(4, @plotColorChartSpectra, expectedWavelengths, expectedSpectra(selectedPatches, :), lineNames, 'expected');
    plotFunWrapper(5, @plotColorChartSpectra, expectedWavelengths, reorderedSpectralValsRaw(selectedPatches, :), lineNames, 'measured-raw');
    plotFunWrapper(6, @plotColorChartSpectra, expectedWavelengths, reorderedSpectralVals(selectedPatches, :), lineNames, 'measured');

    gofs = applyFuncOnRows(reorderedSpectralVals(selectedPatches, :), expectedSpectra(selectedPatches, :), @goodnessOfFit)
    nmses = applyFuncOnRows(reorderedSpectralVals(selectedPatches, :), expectedSpectra(selectedPatches, :), @nmse)


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


end


function [reorderedSpectra, labels] = ReorderSpectra(target, chartColorOrder, spectraColorOrder, wavelengths, spectralWavelengths, isRaw)
if nargin < 6
    isRaw = false;
end
if (isRaw)
    fullReflectance = ones(1, size(target, 2));
else
    load(fullfile(getSetting('matdir'), 'fullreflectance.mat'), 'fullReflectance');
end
[~, idx] = ismember(spectraColorOrder, chartColorOrder);
idx = nonzeros(idx); 
[~, idx2] = ismember(spectralWavelengths', wavelengths);

targetDecim = target(:, idx2);
normConstant = fullReflectance(idx2);

reorderedSpectra = targetDecim(idx, :) ./ normConstant;
labels = spectraColorOrder;
end

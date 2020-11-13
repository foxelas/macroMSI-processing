toRun = false;
if toRun
    prepareData;
end

whites = zeros(length(wavelengths), length(datafiles));
for i = 4:length(datafiles)
    filename = datafiles(i).name;
    [spectralData, imageXYZ, wavelengths] = loadH5Data(filename);

    [colorMasks, chartMask] = getColorCheckerMasks(imageXYZ, filename, true);

    figure(1);
    imagesc(imageXYZ(:, :, 2));
    title('Select ROI for Normalization Spectrum');
    maskWhite = roipoly;
    fullReflectance = getPatchValues(spectralData, maskWhite);
    whites(:, i) = fullReflectance;

    %     croppedLab = xyz2lab(imageXYZ(any(chartMask,2), any(chartMask, 1), :));
    %     [actualLabVals] = getPatchValues(croppedLab, colorMasks);

    croppedSpectral = spectralData(any(chartMask, 2), any(chartMask, 1), :);
    [actualSpectralVals] = getPatchValues(croppedSpectral, colorMasks);
    reorderedSpectralVals = ReorderSpectra(actualSpectralVals, chartColorOrder, spectraColorOrder, wavelengths, expectedWavelengths);
    reorderedSpectralValsRaw = ReorderSpectra(actualSpectralVals, chartColorOrder, spectraColorOrder, wavelengths, expectedWavelengths, true);

    %% Standard (expected) color patch spectra
    plotFunWrapper(4, @plotColorChartSpectra, expectedSpectra, spectraColorOrder, 'expected');

    plotFunWrapper(5, @plotColorChartSpectra, reorderedSpectralValsRaw, spectraColorOrder, 'measured-raw');
    plotFunWrapper(6, @plotColorChartSpectra, reorderedSpectralVals, spectraColorOrder, 'measured');

    gofs = applyFuncOnRows(reorderedSpectralVals, expectedSpectra, @goodnessOfFit)
    nmses = applyFuncOnRows(reorderedSpectralVals, expectedSpectra, @nmse)


    %     %% Get location of color chart
    %     croppedXYZ = imageXYZ(any(chartMask,2), any(chartMask, 1), :);
    %     croppedXYZnorm = croppedXYZ ./ sum(sum(croppedXYZ(:,:,2)));
    %     rgb = xyz2rgb(croppedXYZnorm,'ColorSpace','srgb');
    %
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


function [reorderedSpectra] = ReorderSpectra(target, chartColorOrder, spectraColorOrder, wavelengths, spectralWavelengths, isRaw)
if nargin < 6
    isRaw = false;
end
if (isRaw)
    fullReflectance = ones(1, size(target, 2));
else
    load('fullreflectance.mat', 'fullReflectance');
end
[~, idx] = ismember(spectraColorOrder, chartColorOrder);
[~, idx2] = ismember(spectralWavelengths, wavelengths);
reorderedSpectra = zeros(length(idx), length(idx2));
for i = 1:length(idx)
    for j = 1:length(idx2)
        reorderedSpectra(i, j) = target(idx(i), idx2(j)) ./ fullReflectance(idx2(j));
    end
end

end

disablePrepare = true;
if ~disablePrepare
    prepareData;
end

disable = true; 
allowRoiSelection = false;
setSetting('matdir', 'F:\temp\mspi\matfiles\hsi')

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

[reorderedSpectralVals, lineNames] = reorderSpectra(actualSpectralVals, chartColorOrder, spectraColorOrder, wavelengths, expectedWavelengths);
[reorderedSpectralValsRaw, ~] = reorderSpectra(actualSpectralVals, chartColorOrder, spectraColorOrder, wavelengths, expectedWavelengths);

%% Standard (expected) color patch spectra
plotFunWrapper(4, @plotColorChartSpectra, expectedWavelengths, expectedSpectra(selectedPatches, :), lineNames, {saveFolder, 'expected'});
plotFunWrapper(5, @plotColorChartSpectra, expectedWavelengths, reorderedSpectralVals(selectedPatches, :), lineNames, {saveFolder, 'measured'});

differenceSpectra = expectedSpectra - reorderedSpectralVals;
plotFunWrapper(6, @plotColorChartSpectra, expectedWavelengths, differenceSpectra(selectedPatches, :), lineNames, {saveFolder, 'difference'});

gofs = applyFuncOnRows(reorderedSpectralVals(selectedPatches, :), expectedSpectra(selectedPatches, :), @goodnessOfFit)
nmses = applyFuncOnRows(reorderedSpectralVals(selectedPatches, :), expectedSpectra(selectedPatches, :), @nmse)

%% Standard (expected) color patch spectra after adjustment
[adjustedReorderedSpectralVals, a] = adjustSpectraToWhitePatch(reorderedSpectralVals, spectraColorOrder);

plotFunWrapper(7, @plotColorChartSpectra, expectedWavelengths, adjustedReorderedSpectralVals(selectedPatches, :), lineNames, {saveFolder, 'measured-adjusted'});

differenceSpectraAdjusted = expectedSpectra - adjustedReorderedSpectralVals;
plotFunWrapper(8, @plotColorChartSpectra, expectedWavelengths, differenceSpectraAdjusted(selectedPatches, :), lineNames, {saveFolder, 'difference-adjusted'});

gofs2 = applyFuncOnRows(adjustedReorderedSpectralVals(selectedPatches, :), expectedSpectra(selectedPatches, :), @goodnessOfFit)
nmses2 = applyFuncOnRows(adjustedReorderedSpectralVals(selectedPatches, :), expectedSpectra(selectedPatches, :), @nmse)

lastUpdateDate = date();
matfileSavedir = fullfile(getSetting('savedir'), getSetting('saveFolder'), 'last_run.mat');
fprintf('Saving matfile in %s', matfileSavedir);
save(matfileSavedir, '-v7.3');
disp('Saving matfile finished.');

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




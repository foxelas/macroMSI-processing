function [T, measuredSpectra, adjustedSpectra] = evaluateColorchart(filename, allowRoiSelection, selectedPatches, option)
%EVALUATECOLORCHART returns measured curves end evaluations in comparison
%to expected colorchart curves
%
%   [T, measuredSpectra, adjustedSpectra] = evaluateColorchart(filename,
%   allowRoiSelection, selectedPatches) returns values for
%   evaluation of the colorchart
%
%   Input parameters
%   -filename: HSI file to read
%   -allowRoiSelection: whether manual selection of colorchart ROI is
%   allowed
%   -selectedPatches: indexes of available patches to be used
%   -option: options for normalization when reading the hsi
%
%   Output parameters
%   -T: table including GOF, NMSE, RMSE of the evalueation of measured with
%   expected spectra
%   -measuredSpectra: measured spectra from the HSI
%   -adjustedSpectra: measured spectra from the HSI after adjustment
%

experiment = getSetting('experiment');
patchOrder = getExpectedValues('colorchartOrder', experiment);
[expectedSpectra, patchNames, expectedWavelengths] = getExpectedValues();

if nargin < 2
    allowRoiSelection = false;
end

if nargin < 3
    selectedPatches = 1:length(patchNames);
end

if nargin < 4
    option = [];
end

[spectralData, ~, wavelengths] = loadH5Data(filename, experiment);
[colorMasks, chartMask] = getColorchartMasks(squeeze(spectralData(:, :, 100)), allowRoiSelection, experiment);
actualSpectralVals = readHSI(spectralData, {chartMask, colorMasks}, option);

[reorderedSpectralVals, lineNames] = reorderSpectra(actualSpectralVals, patchOrder, patchNames, wavelengths, expectedWavelengths);
% [reorderedSpectralValsRaw, ~] = reorderSpectra(actualSpectralVals, chartColorOrder, spectraColorOrder, wavelengths, expectedWavelengths);

measuredSpectra = reorderedSpectralVals(selectedPatches, :);
standardSpectra = expectedSpectra(selectedPatches, :);

if ~strcmp(option, 'raw')
    [T, adjustedSpectra] = compareSpectra(standardSpectra, measuredSpectra, lineNames);
end

end
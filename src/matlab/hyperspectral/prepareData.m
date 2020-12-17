
%% main
close all;
%clc;

%% PARAMETERS SET
%Modify userSettings.csv for options
userSettingsFile = '..\..\conf\hsiUserSettings.csv';
readWhite = false;
wavelengths = [380:780]';

% Available dates,  '20201209', '20201111'
dataDate = '20201209';

%% NEXT
setOpt(userSettingsFile);

indir = strcat('F:\temp\mspi\saitama', dataDate, '_test/h5/');
setSetting('datadir', indir);
datafiles = dir(fullfile(indir, '*.h5'));

matfileDir = getSetting('matdir');
m  = matfile(fullfile(matfileDir, 'fullreflectance.mat'), 'Writable', false);

if (~any(strcmp({whos(m).name}, 'maskWhite')) || readWhite)
    i = getIndexValue('whiteReflectance', dataDate);
    %For 99% reflectance
    [spectralData, imageXYZ, wavelengths] = loadH5Data(datafiles(i).name);
    figure(1);
    imagesc(imageXYZ(:, :, 2));
    title('Select ROI for Normalization Spectrum');
    maskWhite = roipoly;
    fullReflectance = getPatchValues(spectralData, maskWhite);
    save('fullreflectance.mat', 'fullReflectance', 'maskWhite');
end
load(fullfile(matfileDir, 'fullreflectance.mat'), 'fullReflectance', 'maskWhite');

setSetting('plotName', fullfile(getSetting('savedir'), 'whitePlot'));
plotFunWrapper(2, @plotSpectra, fullReflectance, wavelengths, '99%-white', 'Reflectance Spectrum of White Balance Sheet');

chartColorOrder = GetOrder();
expectedLabVals = GetExpectedLabValues();
expectedRGBVals = GetExpectedRGBValues();
[expectedSpectra, expectedWavelengths, spectraColorOrder] = GetExpectedSpectra();

function out = GetOrder()
outstruct = delimread(fullfile(getSetting('systemdir'), 'orderOfCurrentChart.txt'), '\t', 'text');
out = outstruct.text;
end

function out = GetExpectedLabValues()
outstruct = delimread(fullfile(getSetting('systemdir'), 'ColorCheckerMicro_Matte_Lab_values.txt'), '\t', 'num');
out = outstruct.num;
end

function out = GetExpectedRGBValues()
outstruct = delimread(fullfile(getSetting('systemdir'), 'ColorCheckerMicro_Matte_RGB_values.txt'), '\t', 'num');
out = outstruct.num;
end

function [spectra, wavelengths, colorNames] = GetExpectedSpectra()
outstruct = delimread(fullfile(getSetting('systemdir'), 'ColorChecker_spectra.txt'), '\t', {'text', 'num'});
colorNames = outstruct.text;
colorNames = colorNames(2:length(colorNames));
wavelengths = outstruct.num(1, :);
spectra = outstruct.num(2:end, :);
end

function searchVal = getIndexValue(searchName, dataDate)
    datasetCharacteristics = delimread( fullfile(getSetting('datasetSettingsDir'), 'dataSetCharacteristics.csv'), ',', {'mixed'}).mixed;
    setId = find(strcmp(datasetCharacteristics(2:end,1), searchName)  &  ([datasetCharacteristics{2:end,2}] == str2num(dataDate))');
    searchVal = datasetCharacteristics(setId + 1, 3);
    searchVal = searchVal{1};
end 

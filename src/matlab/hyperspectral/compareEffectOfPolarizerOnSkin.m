%% Compare effect of polarizing on skin 

close all;

%% Settings 
userSettingsFile = '..\..\conf\hsiUserSettings.csv';

originDir = 'F:\temp\mspi';
dataDate = '20210127';
experiment = 'polarizing_effect_on_tissue';
configuration = 'singleLightClose';
imageFilenames = {'20210127_174840_dry_nofilter_800.h5', ...   
               '20210127_174840_wet_filter_1460.h5', ...    
               '20210127_174840_dry_filter_1460.h5', ...
               '20210127_174840_wet_nofilter_300.h5', ...   
               '20210127_174840_dry_nofilter_1460.h5', ...  
               '20210127_174840_dry_nofilter_300.h5'   };
integrationTimes = [800, 1460, 1460, 300, 1460, 300];

dryInd = [3, 6];
wetInd = [2, 4];

%% Main 
setOpt(userSettingsFile);

indir = fullfile(originDir, '2_saitamaHSI', strcat('saitama', dataDate, '_test'), 'h5');
setSetting('datadir', indir);
matdir = fullfile(originDir, 'matfiles\hsi');
setSetting('matdir', matdir);
setSetting('saveFolder', experiment);
curSaveDir = mkNewDir(getSetting('savedir'), getSetting('saveFolder'));

setSetting('experiment', experiment);
setSetting('configuration', configuration);
setSetting('normByPixel', true);

%% Prepare white
readBlackWhite = false;
if readBlackWhite
    readWhite(dataDate, 300, true, false, 'white', false);
end 

imageFolder = 'skin test';

xPoints = 300:150:800;
yPoints = 300:150:800;

allowRoiSelection = true; 
limits = {[0,0.0035], [0,0.015]};
saveFolderOld = getSetting('saveFolder');
actualSpectralVals = cell(4,1);

%% For dry skin 
baseName = 'dryHand';

for k = 1:length(dryInd) 
    ind = dryInd(k);
    hasFilter = ~contains(imageFilenames{ind}, 'nofilter');
    if hasFilter
        imgName = strcat(baseName, 'Filter');
    else 
        imgName = strcat(baseName, 'NoFilter');
    end
    [filename, integrationTime] = getFilename(configuration, imgName, integrationTimes(ind));
    setSetting('integrationTime', integrationTime);
    
    indir = fullfile(originDir, '2_saitamaHSI', strcat('saitama', dataDate, '_test'), imageFolder, 'h5');
    setSetting('datadir', indir);
    
    [spectralData, imageXYZ, wavelengths] = loadH5Data(filename, experiment);
    baseImage = rescale(spectralData(:,:,200));
    
    setSetting('saveFolder', strcat(saveFolderOld, '\dry', num2str(k)));
    setSetting('plotName',mkNewDir(curSaveDir, imgName));
    plots(1, @plotYFromHSI, baseImage, 'Various points on the image', xPoints, yPoints);
    
    showMultiplePointSpectra(k + 1, xPoints, yPoints, 0, wavelengths', spectralData, limits{k}, imgName);

    [colorMasks, chartMask] = getColorchartMasks(rescale(squeeze(baseImage)), allowRoiSelection, experiment);%(imageXYZ, allowRoiSelection, experiment)
    actualSpectralVals{k} = readHSI(spectralData, {chartMask, colorMasks}, [], hasFilter); 
    lineNames = arrayfun(@num2str, [1:30], 'UniformOutput', false);
    plots(1, @plotColorChartSpectra, actualSpectralVals{k}, lineNames, 'measured-points');
end
setSetting('saveFolder', saveFolderOld); 

%% For wet skin
close all;
baseName = 'wetHand';

for k = 1:length(wetInd) 
    ind = wetInd(k);
    hasFilter = ~contains(imageFilenames{ind}, 'nofilter');
    if hasFilter
        imgName = strcat(baseName, 'Filter');
    else 
        imgName = strcat(baseName, 'NoFilter');
    end
    [filename, integrationTime] = getFilename(configuration, imgName, integrationTimes(ind));
    setSetting('integrationTime', integrationTime);
    
    indir = fullfile(originDir, '2_saitamaHSI', strcat('saitama', dataDate, '_test'), imageFolder, 'h5');
    setSetting('datadir', indir);
    
    [spectralData, imageXYZ, wavelengths] = loadH5Data(filename, experiment);
    baseImage = rescale(spectralData(:,:,200));
    
    setSetting('saveFolder', strcat(saveFolderOld, '\wet', num2str(k)));
    setSetting('plotName',mkNewDir(curSaveDir, imgName));
    plots(1, @plotYFromHSI, baseImage, 'Various points on the image', xPoints, yPoints);
    
    showMultiplePointSpectra(k + 1, xPoints, yPoints, 0, wavelengths', spectralData, limits{k}, imgName);

    [colorMasks, chartMask] = getColorchartMasks(rescale(squeeze(baseImage)), allowRoiSelection, experiment);%(imageXYZ, allowRoiSelection, experiment)
    actualSpectralVals{k + 2}  = readHSI(spectralData, {chartMask, colorMasks}, [], hasFilter); 
    lineNames = arrayfun(@num2str, [1:30], 'UniformOutput', false);
    plots(1, @plotColorChartSpectra, actualSpectralVals{k}, lineNames, 'measured-points');
end
setSetting('saveFolder', saveFolderOld); 

savedir = fullfile(getSetting('savedir'), getSetting('saveFolder'), 'normalizedSpectra.mat');
save( savedir, 'actualSpectralVals');
fprintf('Save actual spectral values in %s\n\n', savedir);

%% Show images with error 

setSetting('saveImages', false);
close all; 
%800
[filename, ~] = getFilename(configuration, 'dryHandNoFilter', 800);

indir = fullfile(originDir, '2_saitamaHSI', strcat('saitama', dataDate, '_test'), imageFolder, 'h5');
setSetting('datadir', indir);

[spectralData, ~, wavelengths] = loadH5Data(filename, experiment);
baseImage = rescale(spectralData(:,:,200));
figure(1)
imagesc(baseImage);
showMultiplePointSpectra(2, xPoints, yPoints, 0, wavelengths', spectralData, limits, imgName);
ylim([0,0.01]);

%1460
[filename, integrationTime] = getFilename(configuration, 'dryHandNoFilter', 1460);

[spectralData, ~, wavelengths] = loadH5Data(filename, experiment);
baseImage = rescale(spectralData(:,:,200));
figure(3)
imagesc(baseImage);

showMultiplePointSpectra(4, xPoints, yPoints, 0, wavelengths', spectralData, limits, imgName);
ylim([0,0.01]);
    
    
    
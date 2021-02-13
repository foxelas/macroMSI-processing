%% Compare image averaging 

close all;

%% Settings 
userSettingsFile = '..\..\conf\hsiUserSettings.csv';

originDir = 'F:\temp\mspi';
dataDate = '20210127';
experiment = 'capture_average_comparison';
configuration = 'singleLightClose';
integrationTime = 1460;
toFuseNames = {'20210127_122055_top_left_1.h5', '20210127_122055_top_left_2.h5', '20210127_122055_top_left_2b.h5'};
baseName = 'colorchartTopLeft';

%% Main 
setOpt(userSettingsFile);

indir = fullfile(originDir, '2_saitamaHSI', strcat('saitama', dataDate, '_test'), 'h5');
setSetting('datadir', indir);
matdir = fullfile(originDir, 'matfiles\hsi');
setSetting('matdir', matdir);

n = length(toFuseNames);

setSetting('integrationTime', integrationTime);
setSetting('experiment', experiment);
setSetting('configuration', configuration);
setSetting('normByPixel', true);
%% Prepare white
readBlackWhite = false;
if readBlackWhite
    readWhite(dataDate,integrationTime, true, false, 'white');
end 

%% Select various points on the image and visualize spectra 
if false 
setSetting('saveFolder', strcat(experiment, '\singlePixelMask')); %'\squarePatchMask'
curSaveDir = mkNewDir(getSetting('savedir'), getSetting('saveFolder'));

spectVals =  zeros(n, 6, 5, 401);
for k = 1:n
    
    imgName = strcat(baseName, num2str(k));
    filename = getFilename(configuration, imgName);
    indir = fullfile(originDir, '2_saitamaHSI', strcat('saitama', dataDate, '_test'),'image fusion', 'h5');
    setSetting('datadir', indir);
    
    [spectralData, imageXYZ, wavelengths] = loadH5Data(filename, experiment);
    
    xPoints = 50:200:1088;
    yPoints = 50:200:982;
    
%     %% plot Y image only
%     setSetting('plotName',mkNewDir(curSaveDir, imgName));
%     plots(1, @plotYFromHSI, imageXYZ, 'Luminance Y image');
    
    %% plot Y image with various points     
    setSetting('plotName',mkNewDir(curSaveDir, imgName));
    plots(2, @plotYFromHSI, imageXYZ, 'Various points on the image', xPoints, yPoints);
    
    limits = [0,0.008];
    spectVals(k, :, :, :) = reshape(showMultiplePointSpectra(k, xPoints, yPoints, 0, wavelengths', spectralData, limits), [1, numel(xPoints), numel(yPoints), numel(wavelengths)]);

end

spectValsAverage = squeeze(mean(spectVals, 1));
showMultiplePointSpectra(4, 1:numel(xPoints), 1:numel(yPoints), 0, wavelengths, spectValsAverage, limits, strcat(imgName, '3captureAverage'));

% spectVals(4, :, :, :) = spectValsAverage;
% spectVals = reshape(spectVals, [size(spectVals, 1), numel(xPoints) * numel(yPoints), length(wavelengths)]);
showMultiplePointSpectra(5, xPoints, yPoints, 0, wavelengths, spectVals, limits, strcat(imgName, 'allResults'));
end 

%% Compare curves 
close all; 
tables = cell(n+1, 1);
measuredSpectra = cell(n+1,1);
adjustedSpectra = cell(n+1,1);
allowRoiSelection = true;
avgMeasured = zeros(24, 36);

for k = 1:n
    
    imgName = strcat(baseName, num2str(k));
    filename = getFilename(configuration, imgName);
    indir = fullfile(originDir, '2_saitamaHSI', strcat('saitama', dataDate, '_test'),'image fusion', 'h5');
    setSetting('datadir', indir);
    
    [tables{k}, measuredSpectra{k}, adjustedSpectra{k}] = evaluateColorchart(filename, experiment, allowRoiSelection);
    avgMeasured = avgMeasured + measuredSpectra{k};
end
measuredSpectra{n+1} = avgMeasured ./ n; 
expected = getExpectedValues();
tables{n+1} = compareSpectra(expected, measuredSpectra{n+1}, [tables{k}.Patch]);

for i = 1:(n+1)
tables{i,1}{25,1} = {'Average'};
tables{i,1}{26,1} = {'StandardDeviation'};
    for j = 2:7
    tables{i,1}{25,j} = mean(tables{i,1}{:,j});
    tables{i,1}{26,j} = std(tables{i,1}{:,j});
    end
end 

for j = 1:24
    if tables{1,1}.AdjGoF(j) < 0.8
        fprintf('**Low result for patch %s\n', tables{1,1}.Patch{j});
    end
end
    
savedir = fullfile(getSetting('savedir'), getSetting('saveFolder'), 'resultTable.mat');
save(savedir, 'tables', 'measuredSpectra', 'adjustedSpectra');
fprintf('Saved values in %s,\n', savedir);

for i = 1:(n+1)
    writetable(tables{i,1}, fullfile(getSetting('savedir'), getSetting('saveFolder'), 'resultTable.xlsx'),'Sheet',i);
end

%% For measured 
close all; 
multipleSpectra = zeros(size(measuredSpectra{1},1), size(measuredSpectra{1},2), n);
for i = 1:n
    multipleSpectra(:,:,i) = measuredSpectra{i};
end 
[~, patchNames] = getExpectedValues();
plots(1, @plotColorChartSpectra, multipleSpectra, patchNames, 'measured_imageAverage');

%% For adjusted 
close all; 
multipleSpectra = zeros(size(adjustedSpectra{1},1), size(adjustedSpectra{1},2), n);
for i = 1:n
    multipleSpectra(:,:,i) = adjustedSpectra{i};
end 
plots(1, @plotColorChartSpectra, multipleSpectra, patchNames,'measured-adjusted_imageAverage');

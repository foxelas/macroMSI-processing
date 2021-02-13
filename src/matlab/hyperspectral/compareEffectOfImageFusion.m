%% Compare image fusion 
close all;

%% Settings 
userSettingsFile = '..\..\conf\hsiUserSettings.csv';

originDir = 'F:\temp\mspi';
dataDate = '20210127';
experiment = 'fusion_comparison';
configuration = 'singleLightClose';
baseName = 'colorchartTopRight'; %'colorchartTopLeft' 'colorchartBottomRight''colorchartBottomLeft'

%% Main 
setOpt(userSettingsFile);

indir = fullfile(originDir, '2_saitamaHSI', strcat('saitama', dataDate, '_test'), 'h5');
setSetting('datadir', indir);
matdir = fullfile(originDir, 'matfiles\hsi');
setSetting('matdir', matdir);
n = length(8);

setSetting('experiment', experiment);
setSetting('configuration', configuration);
setSetting('normByPixel', true);

%% Prepare white
readBlackWhite = false;
if readBlackWhite
    for i = [ 500, 800, 1460] %1800,
        integrationTime = i;
        setSetting('integrationTime', integrationTime);
        readWhite(dataDate,integrationTime, true, false, experiment, configuration, 'white');
    end
end 

%% Select various points on the image and visualize spectra 
if false 
spectVals =  zeros(n, 6, 5, 401);
for k = 1:n
    
    imgName = strcat(baseName, num2str(k));
    filename = getFilename(configuration, imgName);
    indir = fullfile(originDir, '2_saitamaHSI', strcat('saitama', dataDate, '_test'),'image fusion', 'h5');
    setSetting('datadir', indir);
    
    [spectralData, imageXYZ, wavelengths] = loadH5Data(filename, experiment);
    
    xPoints = 50:200:1088;
    yPoints = 50:200:982;
    
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

ind = 0;  
for k = 1:9
    if k ~= 3
        ind = ind + 1;
        imgName = strcat(baseName, num2str(k));
        [filename, integrationTime]= getFilename(configuration, imgName);
        indir = fullfile(originDir, '2_saitamaHSI', strcat('saitama', dataDate, '_test'),'image fusion', 'h5');
        setSetting('datadir', indir);
        setSetting('saveFolder', strcat(experiment, '\singlePixelMask\', baseName, '\', num2str(k)));
         
        [tables{ind}, measuredSpectra{ind}, adjustedSpectra{ind}] = evaluateColorchart(filename, experiment, allowRoiSelection);
    end 
end

fusedSpectra = cell(6,1);

%% Single image at 1460ms
avgMeasured1 = measuredSpectra{1};
fusedSpectra{1} = avgMeasured1;

%% Fuse 2 images at 1460ms 
avgMeasured2 = (measuredSpectra{1} + measuredSpectra{2}) /2;
fusedSpectra{2} = avgMeasured2;

%% Single image at 1800ms
avgMeasured3 = measuredSpectra{3};
fusedSpectra{3} = avgMeasured3;

%% Fuse 2 images at 1460 (380-540nm) and 1800 (541-780nm)
[~,~,ia] = intersect(expectedWavelengths, [380:540]);
% [~,~,ib] = intersect(expectedWavelengths, [541:780]);

avgMeasured4 = zeros(24, 36);
avgMeasured4(:, 1:length(ia)) = measuredSpectra{4};
avgMeasured4(:, length(ia)+1:end) = measuredSpectra{5};
fusedSpectra{4} = avgMeasured4;

%% Fusion image at 1460 and 1800 
avgMeasured5 = zeros(24, 36);
avgMeasured5(:, 1:length(ia)) = measuredSpectra{1}(:, 1:length(ia));
avgMeasured5(:, length(ia)+1:end) = (measuredSpectra{1}(:, length(ia)+1:end) + measuredSpectra{5}) /2;
fusedSpectra{5} = avgMeasured5;

%% Fusion image at 2x500 ms (380-780nm) and 1x800ms (541-780nm)
avgMeasured6 = zeros(24, 36);
avgMeasured6 = (measuredSpectra{6} + measuredSpectra{7}) / 2;
avgMeasured6(:, length(ia)+1:end) = (avgMeasured6(:, length(ia)+1:end) + measuredSpectra{5})/ 2;
fusedSpectra{6} = avgMeasured6;

n = length(fusedSpectra);
fusedSpectraTables = cell(n, 1);
[expected, patchNames] = getExpectedValues();

for i = 1:n
    setSetting('saveFolder', strcat(experiment, '\singlePixelMask\', baseName, '\fusedSpectra', num2str(i)));
    fusedSpectraTables{i} = compareSpectra(expected, ...
        fusedSpectra{i}, [tables{1}.Patch]);
end

for i = 1:n
fusedSpectraTables{i,1}{25,1} = {'Average'};
fusedSpectraTables{i,1}{26,1} = {'StandardDeviation'};
    for j = 2:7
    fusedSpectraTables{i,1}{25,j} = mean(fusedSpectraTables{i,1}{:,j});
    fusedSpectraTables{i,1}{26,j} = std(fusedSpectraTables{i,1}{:,j});
    end
end 

for j = 1:24
    if fusedSpectraTables{1,1}.AdjGoF(j) < 0.8
        fprintf('**Low result for patch %s\n', fusedSpectraTables{1,1}.Patch{j});
    end
end


setSetting('saveFolder', strcat(experiment, '\singlePixelMask\', baseName));
savedir = fullfile(getSetting('savedir'), getSetting('saveFolder'), 'resultTable.mat');
save(savedir, 'tables', 'measuredSpectra', 'adjustedSpectra', 'fusedSpectra', 'fusedSpectraTables');
fprintf('Saved values in %s,\n', savedir);

for i = 1:n
    writetable(fusedSpectraTables{i,1}, fullfile(getSetting('savedir'), getSetting('saveFolder'), 'resultTable.xlsx'),'Sheet',i);
end

%% For measured 
close all; 
multipleSpectra = zeros(size(fusedSpectra{1},1), size(fusedSpectra{1},2), n);
for i = 1:n
    multipleSpectra(:,:,i) = fusedSpectra{i};
end 
plots(1, @plotColorChartSpectra, multipleSpectra, patchNames, 'measured_fused');

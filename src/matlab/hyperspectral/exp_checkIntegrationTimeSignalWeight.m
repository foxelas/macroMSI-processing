%% Set the colorchart at different positions on the capture base and 
% compare with the measured values from source closer to the light source
% and at different integration times 

%% Setup 
startRun;

%% Colorchart with different normalizations and positions
experiment = 'testIntegrationTimeSignalWeight';
setSetting('experiment', experiment);
dataDate = '20210127';
configuration = 'singleLightClose';
normalization = 'raw';
initialization;

integrationTimes = [1460, 1800, 500, 800];
z = numel(integrationTimes);

for j = 1:z
    integrationTime = integrationTimes(j);
    setSetting('integrationTime', integrationTime);
    setSetting('normalization', normalization);
    setSetting('saveFolder', fullfile(experiment, num2str(integrationTime)));
    target = 'white_filter';
    readHSIData('whiteReflectance', target, experiment);
end

xPoints = [50, 1350, 700,  50, 1350] + 2;
yPoints = [50, 50, 500, 975, 975] + 2;

pointMasks = zeros(1376, 1024,5);
for i = 1:5
    pointMasks(xPoints(i), yPoints(i),i) = 1; 
end 

raws = cell(z, 1);
spects = cell(z,1);
for j = 1:z
    integrationTime = integrationTimes(j);
    setSetting('integrationTime', integrationTime);
    setSetting('normalization', normalization);
    setSetting('saveFolder', fullfile(experiment, num2str(integrationTime)));
    target = 'white_filter';
    
    fileConditions = getFileConditions('whiteReflectance', target);
    filename = getFilename(fileConditions{:});
    [raw, ~, wavelengths] = loadH5Data(filename);
    raws{j} = raw;
    spects{j} = getSpectrumCurves(raw, pointMasks);
end

for i = 1:5
    pointVersions = zeros(z, 401);
    lineNames = cell(z, 1);
    for j = 1:z 
        pointVersions(j,:) = spects{j,1}(i,:);
        lineNames{j} = num2str(integrationTimes(j));
    end 
    setSetting('saveFolder', fullfile(experiment, strcat('point_', num2str(xPoints(i)), '_', num2str(yPoints(i)))));
    plots(1, @plotColorChartSpectra, pointVersions, lineNames, 'measured-raw', [0, 0.003]);
end 

baseImage = getDisplayImage(raw, 'rgb');
plotName = fullfile(getSetting('savedir'), getSetting('experiment'), strcat(target, '-points', '.png'));
setSetting('plotName', plotName);
plots(2, @plotPointsOnImage, baseImage, yPoints, xPoints, false);

idxs = {'1', '3', '6'};
z  = 3;
selectedPatchIndex = [2, 13, 14, 15]; %19
allowRoiSelection = true;
patchSpects = cell(z,1);
for j = 1:z
    integrationTime = integrationTimes(j);
    setSetting('integrationTime', integrationTime);
    setSetting('normalization', 'raw');
    target = strcat('colorchartBottomLeft', idxs{j});
    readHSIData('colorchart', target, experiment);
    setSetting('saveFolder', fullfile(experiment, strcat(target, '_', num2str(integrationTime))));

    [~, patchSpects{j}, ~, ~] = evaluateColorchart(target, allowRoiSelection, selectedPatchIndex, expectedSpectra); 
end

lineColor = {'lightSkin', 'blue', 'green', 'red'};
k = numel(selectedPatchIndex);
for i = 1:k
    pointVersions = zeros(z, 36);
    lineNames = cell(z, 1);
    for j = 1:z 
        pointVersions(j,:) = patchSpects{j,1}(i,:);
        lineNames{j} = num2str(integrationTimes(j));
    end 
    setSetting('saveFolder', fullfile(experiment, strcat('point_', lineColor{i})));
    plots(1, @plotColorChartSpectra, pointVersions, lineNames, 'measured-raw', [0, 0.003]);
end 


k = numel(selectedPatchIndex);
z  = 3;
rmses = zeros(k, z);
for i = 1:k
    pointVersions = zeros(z, 36);
    lineNames = cell(z, 1);
    for j = 1:z 
        pointVersions(j,:) = patchSpects{j,1}(i,:);
        lineNames{j} = num2str(integrationTimes(j));
    end 
    
    rmses(i,1) = rmse(pointVersions(1,:),pointVersions(2,:));
    rmses(i,2) = rmse(pointVersions(3,:), pointVersions(2,:));
    rmses(i,3) = rmse(pointVersions(3,:), pointVersions(1,:));
end 

Tpatch = array2table(rmses,'VariableNames',{'1460vs1800', '500vs1800', '500vs1460'});
Tpatch.Patch = deal(lineColor');

names = cell(5,1);
z = 4;
rmses = zeros(k, 5);
for i = 1:5
    pointVersions = zeros(z, 401);
    lineNames = cell(z, 1);
    for j = 1:z 
        pointVersions(j,:) = spects{j,1}(i,:);
        lineNames{j} = num2str(integrationTimes(j));
    end 
    
    rmses(i,1) = rmse(pointVersions(1,:),pointVersions(2,:));
    rmses(i,2) = rmse(pointVersions(3,:), pointVersions(2,:));
    rmses(i,3) = rmse(pointVersions(3,:), pointVersions(1,:));
    rmses(i,4) = rmse(pointVersions(4,:), pointVersions(1,:));
    rmses(i,5) = rmse(pointVersions(3,:), pointVersions(4,:));

    names{i} = strcat('point_', num2str(xPoints(i)), '_', num2str(yPoints(i)));
end 


Twhite = array2table(rmses,'VariableNames',{'1460vs1800', '500vs1800', '500vs1460', '800vs1460', '500vs800'});
Twhite.Patch = deal(names);

Twhite

Tpatch

v = fullfile(getSetting('savedir'), getSetting('experiment'));
save(strcat(v, 'lastrun2.mat'), '-v7.3');
endRun;

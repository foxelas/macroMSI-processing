%% Set the colorchart at different positions on the capture base and 
% compare with Babel expected values 

%% Setup 
startRun;

%% Colorchart with different normalizations and positions

% version = '2'; 
version = '3';
if strcmp(version, '2')
    dataDate = '20210317';
    integrationTime = 1360;
    colorPatchOrder = 'bluishGreenRight';

elseif strcmp(version, '3')
    dataDate ='20210406';
    integrationTime = 618;
    colorPatchOrder = 'redRight';
end

experiment = strcat('testCalibrationPositions', version);
setSetting('experiment', experiment);
configuration = 'singleLightClose';
normalization = 'byPixel';
initialization;

setSetting('isRotated', false);
allowRoiSelection = true;
positions = {'BottomLeft', 'TopLeft', 'TopRight','BottomRight', 'Middle'};
n = numel(positions);
for i = 1:n
    setSetting('normalization', 'byPixel');
    target = strcat('colorchart', positions{i});
    readHSIData('colorchart', target, experiment);
end 

% warning('Running for patches: light skin, red, green, blue');
% selectedPatchIndex = [2, 13, 14, 15];
selectedPatchIndex = 1:24;

for i = 1:n
    target = strcat('colorchart', positions{i});
    setSetting('saveFolder', fullfile(experiment, 'raw',  positions{i}));
    setSetting('normalization', 'raw');
    evaluateColorchart(target, allowRoiSelection, selectedPatchIndex); 
end 

normalizations = { 'bandmax', 'uniSpectrum', 'byPixel'};
m = numel(normalizations);
tables = cell(m,1);
measuredSpectra = cell(m,1);
adjustedSpectra = cell(m,1);
alphas = cell(m,1);
confs = cell(m,2);

k = 0;
for i = 1:n
    for j = 1:m
        k = k + 1;
        normalization = normalizations{j};
        position =  positions{i};
        target = strcat('colorchart', position);
        setSetting('saveFolder', fullfile(experiment, normalization, position));
        setSetting('normalization', normalization);
        confs(k, 1:2) = deal({configuration, strcat(normalization, '_', position)});
  
        [tables{k}, measuredSpectra{k}, adjustedSpectra{k}, alphas{k}] = evaluateColorchart(target, allowRoiSelection, selectedPatchIndex); 
    end 
end 

[tables] = exportSimilarityTables(tables, measuredSpectra, adjustedSpectra, alphas, confs);

endRun;
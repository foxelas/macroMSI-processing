%% Set the colorchart at different positions on the capture base and 
% compare with Babel expected values 

%% Setup 
startRun;

%% Colorchart with different normalizations and positions
experiment = 'testCalibrationPositions';
setSetting('experiment', experiment);
dataDate = '20210317'; %'20210406';
integrationTime = 1360; % 618;
configuration = 'singleLightClose';
normalization = 'byPixel';
initialization;

setSetting('colorPatchOrder',  'bluishGreenRight'); %'redRight'
setSetting('isRotated', false);
allowRoiSelection = true;
positions = {'BottomLeft', 'TopLeft', 'TopRight','BottomRight', 'Middle'};
n = numel(positions);
for i = 1:n
    setSetting('normalization', 'byPixel');
    target = strcat('colorchart', positions{i});
    readHSIData('colorchart', target, experiment);
end 

warning('Running for patches: light skin, red, green, blue');
selectedPatchIndex = [2, 13, 14, 15];

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
%% Set the colorchart at different positions on the capture base and 
% compare with the measured values from source closer to the light source 

%% Setup 
startRun;

%% Colorchart with different normalizations and positions
experiment = 'testCalibrationPositionsRelative';
setSetting('experiment', experiment);
dataDate = '20210317';
integrationTime = 1360;
configuration = 'singleLightClose';
normalization = 'byPixel';
initialization;

setSetting('colorPatchOrder', 'bluishGreenRight');
setSetting('isRotated', false);
allowRoiSelection = true;
selectedPatchIndex = [2, 13, 14, 15]; %19
positions = {'BottomLeft', 'TopLeft', 'TopRight','BottomRight', 'Middle'};
n = numel(positions);
normalizations = {'bandmax', 'byPixel'};

m = numel(normalizations);
tables = cell(m,1);
measuredSpectra = cell(m,1);
adjustedSpectra = cell(m,1);
alphas = cell(m,1);
confs = cell(m,2);

for i = 1:n 
    position =  positions{i};
    target = strcat('colorchart', position);
    readHSIData('colorchart', target, experiment);
    setSetting('normalization', 'raw');
    setSetting('saveFolder', fullfile(experiment, 'raw', position));
    evaluateColorchart(target, allowRoiSelection, selectedPatchIndex); 
end 

k = 0;
for j = 1:m
    normalization = normalizations{j};
    setSetting('normalization', normalization);
    for i = 1:n
        position =  positions{i};
        target = strcat('colorchart', position);
        setSetting('saveFolder', fullfile(experiment, normalization, position));
        
        if i ~= 1
            k = k + 1;
            confs(k, 1:2) = deal({configuration, strcat(normalization, '_', position)});
            [tables{k}, measuredSpectra{k}, adjustedSpectra{k}, alphas{k}] = evaluateColorchart(target, allowRoiSelection, selectedPatchIndex, expectedSpectra); 
        else
            [~, expectedSpectra, ~, ~] = evaluateColorchart(target, allowRoiSelection, selectedPatchIndex); 
        end
    end 

end 

[tables] = exportSimilarityTables(tables, measuredSpectra, adjustedSpectra, alphas, confs);

endRun;
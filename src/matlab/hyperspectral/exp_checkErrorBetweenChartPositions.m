%% Set the colorchart at different positions on the capture base and 
% compare with the measured values from source closer to the light source 
%   To run after scipt: exp_checkErrorAtDifferentChartPositions.m 

%% Setup 
startRun;

%% Colorchart with different normalizations and positions
experiment = 'testCalibrationPositionsRelative';
setSetting('experiment', experiment);
dataDate = '20210317'; %'20210406';
integrationTime = 1360; % 618;
configuration = 'singleLightClose';
normalization = 'byPixel';
initialization;

setSetting('colorPatchOrder',  'bluishGreenRight'); %'redRight'
setSetting('isRotated', false);
allowRoiSelection = true;
selectedPatchIndex = [2, 13, 14, 15]; %19
positions = {'Middle', 'BottomLeft', 'TopLeft', 'TopRight','BottomRight'};
n = numel(positions);
normalizations = {'bandmax', 'byPixel'};

m = numel(normalizations);
tables = cell(m,1);
measuredSpectra = cell(m,1);
adjustedSpectra = cell(m,1);
alphas = cell(m,1);
confs = cell(m,2);

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
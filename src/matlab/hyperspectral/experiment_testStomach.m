%% Setup 
startRun;

%% Read tissue 1  
version = '1';
experiment = strcat('testStomach', version);
dataDate = '20210308';
integrationTime = 1460;
configuration = 'singleLightClose';
initialization;

% normalizations = {'raw', 'bandmax', 'uniSpectrum', 'byPixel'};
% m = numel(normalizations);
% target = strcat('stomachTissue', version);
% xPoints = [250, 320, 400];
% yPoints = [250, 350, 400];
% % readHSIData('tissue', target, experiment);
% setSetting('isRotated', false);
% for k = 1:m
%     normalization = normalizations{k};
%     setSetting('saveFolder', fullfile(experiment, normalization));
%     setSetting('normalization', normalization);    
%     fileConditions = getFileConditions('tissue', target);
%     [measured, curveNames] = getRepresentativePoints(fileConditions, xPoints, yPoints);
% end

%% Colorchart middle with different normalizations
setSetting('experiment', experiment);
setSetting('isRotated', false);
readHSIData('colorchart', target, experiment);
allowRoiSelection = true;
setSetting('colorPatchOrder', 'bluishGreenRight');

setSetting('saveFolder', fullfile(experiment, 'raw'));
setSetting('normalization', 'raw');
evaluateColorchart(target, allowRoiSelection); 
    
normalizations = { 'bandmax', 'uniSpectrum', 'byPixel'};
m = numel(normalizations);
m = numel(normalizations);
tables = cell(m,1);
measuredSpectra = cell(m,1);
adjustedSpectra = cell(m,1);
alphas = cell(m,1);
confs = cell(m,2);
for k = 1:m
    normalization = normalizations{k};
    setSetting('saveFolder', fullfile(experiment, normalization));
    setSetting('normalization', normalization);
    confs(k, 1:2) = deal({configuration, normalization});
    
    [tables{k}, measuredSpectra{k}, adjustedSpectra{k}, alphas{k}] = evaluateColorchart(target, allowRoiSelection); 
end 

[tables] = exportSimilarityTables(tables, measuredSpectra, adjustedSpectra, alphas, confs);

endRun;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read tissue 2
version = '2';
experiment = strcat('testStomach', version);
dataDate = '20210317';
integrationTime = 1360;
configuration = 'singleLightClose';
normalization = 'byPixel';

initialization;

normalizations = {'raw', 'bandmax', 'uniSpectrum', 'byPixel'};
m = numel(normalizations);
target = strcat('stomachTissue', version);
xPoints = [150, 200, 250];
yPoints = [250, 325, 400];
readHSIData('tissue', target, experiment);
setSetting('isRotated', false);
for k = 1:m
    normalization = normalizations{k};
    setSetting('saveFolder', fullfile(experiment, normalization));
    setSetting('normalization', normalization);    
    fileConditions = getFileConditions('tissue', target);
    [measured, curveNames] = getRepresentativePoints(fileConditions, xPoints, yPoints);
end 

%% Colorchart middle with different normalizations
experiment = 'testCalibration';
setSetting('experiment', experiment);
setSetting('isRotated', false);
target = strcat('colorchart', 'Middle');
readHSIData('colorchart', target, experiment);
allowRoiSelection = true;
setSetting('colorPatchOrder', 'bluishGreenRight');

setSetting('saveFolder', fullfile(experiment, 'raw'));
setSetting('normalization', 'raw');
evaluateColorchart(target, allowRoiSelection); 
    
normalizations = { 'bandmax', 'uniSpectrum', 'byPixel'};
m = numel(normalizations);
m = numel(normalizations);
tables = cell(m,1);
measuredSpectra = cell(m,1);
adjustedSpectra = cell(m,1);
alphas = cell(m,1);
confs = cell(m,2);
for k = 1:m
    normalization = normalizations{k};
    setSetting('saveFolder', fullfile(experiment, normalization));
    setSetting('normalization', normalization);
    confs(k, 1:2) = deal({configuration, normalization});
    
    [tables{k}, measuredSpectra{k}, adjustedSpectra{k}, alphas{k}] = evaluateColorchart(target, allowRoiSelection); 
end 

[tables] = exportSimilarityTables(tables, measuredSpectra, adjustedSpectra, alphas, confs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Colorchart with different normalizations and positions
experiment = 'testCalibrationPositions';
setSetting('experiment', experiment);

setSetting('colorPatchOrder', 'bluishGreenRight');
setSetting('isRotated', false);
allowRoiSelection = true;
positions = {'BottomLeft', 'TopLeft', 'TopRight','BottomRight', 'Middle'};
n = numel(positions);
for i = 1:n
    target = strcat('colorchart', positions{i});
    readHSIData('colorchart', target, experiment);

    setSetting('saveFolder', fullfile(experiment, 'raw',  positions{i}));
    setSetting('normalization', 'raw');
    evaluateColorchart(target, allowRoiSelection); 
end 

normalizations = { 'bandmax', 'uniSpectrum', 'byPixel'};
m = numel(normalizations);
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
  
        [tables{k}, measuredSpectra{k}, adjustedSpectra{k}, alphas{k}] = evaluateColorchart(target, allowRoiSelection); 
    end 
end 

[tables] = exportSimilarityTables(tables, measuredSpectra, adjustedSpectra, alphas, confs);

endRun;



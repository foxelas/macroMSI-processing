%% Set stomach tissue and colorchart at middle for March 17th 
% compare with the measured values from Babel and for different
% normalization

%% Setup 
startRun;

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
experiment = 'testCalibration2';
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

endRun;
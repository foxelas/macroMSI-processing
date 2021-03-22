%% Setup 
startRun;

version = '2';
experiment = strcat('testStomach', version);
if strcmp(version, '1')
    dataDate = '20210308';
    integrationTime = 1460;
else 
    dataDate = '20210317';
    integrationTime = 1360;
end 

configuration = 'singleLightClose';
normalization = 'byPixel';

initialization;

setSetting('saveFolder', experiment);

%% Read white and black
readWhite(dataDate, integrationTime, experiment, configuration, []);

%% Read colorchart 
normalizations = {'bandmax', 'uniSpectrum', 'byPixel'};
m = numel(normalizations);
tables = cell(m,1);
measuredSpectra = cell(m,1);
adjustedSpectra = cell(m,1);
alphas = cell(m,1);
allowRoiSelection = false;
confs = cell(m,2);
target = strcat('stomachTissue', version);

setSetting('isRotated', false);
for k = 1:m
    normalization = normalizations{k};
    setSetting('saveFolder', fullfile(experiment, normalization));
    setSetting('normalization', normalization);
    confs(k, 1:2) = deal({configuration, normalization});
    fileConditions = getFileConditions('colorchart', target);
    [tables{k}, measuredSpectra{k}, adjustedSpectra{k}, alphas{k}] = evaluateColorchart(fileConditions, allowRoiSelection); 

    %% Read sample 
    xPoints = [250, 320, 400];
    yPoints = [250, 350, 400];
    fileConditions = getFileConditions('tissue', target);
    [measured, curveNames] = getRepresentativePoints(fileConditions, xPoints, yPoints);
end 

[tables] = exportSimilarityTables(tables, measuredSpectra, adjustedSpectra, alphas, confs);
endRun;



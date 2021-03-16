%% Setup 
startRun;

experiment = 'testStomach';
dataDate = '20210308';
configuration = 'singleLightClose';
integrationTime = 1460;
normalization = 'byPixel';

initialization;

setSetting('saveFolder', experiment);

%% Read white and black
% readWhite(dataDate, integrationTime, experiment, configuration, []);

%% Read colorchart 
normalizations = {'bandmax', 'uniSpectrum', 'byPixel'};
m = numel(normalizations);
tables = cell(m,1);
measuredSpectra = cell(m,1);
adjustedSpectra = cell(m,1);
alphas = cell(m,1);
allowRoiSelection = false;
confs = cell(m,2);


setSetting('isRotated', false);
for k = 1:m
    normalization = normalizations{k};
    setSetting('saveFolder', fullfile(experiment, normalization));
    setSetting('normalization', normalization);
    confs(k, 1:2) = deal({configuration, normalization});
    [tables{k}, measuredSpectra{k}, adjustedSpectra{k}, alphas{k}] = evaluateColorchart('colorchart', allowRoiSelection); 

    %% Read sample 
    xPoints = [250, 320, 400];
    yPoints = [250, 350, 400];
    [measured, curveNames] = getRepresentativePoints('stomachSample', xPoints, yPoints);
end 

[tables] = exportSimilarityTables(tables, measuredSpectra, adjustedSpectra, alphas, confs);
endRun;

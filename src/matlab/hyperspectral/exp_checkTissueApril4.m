%% Set stomach tissue and colorchart at middle 
% compare with the measured values from Babel and for different
% normalization

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% for April 6th %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup 
startRun;

%% Read tissue 2
version = '3';
experiment = strcat('testStomach', version);
dataDate = '20210406';

targets = { strcat('stomachTissue', version), strcat('stomachTissue', version, '_int_auto_gain_4')};
integrationTimes = [618, 160];
types = {'tissue', 'tissue'};

targetInfo.Targets = targets;
targetInfo.IntegartionTimes = integrationTimes;
targetInfo.Types = types; 
importImagesToMat(dataDate, experiment, targetInfo);


normalizations = {'raw', 'bandmax', 'uniSpectrum', 'byPixel'};
m = numel(normalizations);
xPoints = [500, 600, 700];
yPoints = [500, 600, 700];

for k = 1:m
    for i = 1:2
        target = targets{i};
        setSetting('integrationTime', integrationTimes(i))
        normalization = normalizations{k};
        setSetting('saveFolder', fullfile(experiment, normalization));
        setSetting('normalization', normalization);    
        fileConditions = getFileConditions('tissue', target);
        [measured, curveNames] = getRepresentativePoints(fileConditions, xPoints, yPoints);
    end
end 

target = 'fullscreen';
setSetting('normalization', 'raw');
integrationTimes = [1380, 800, 2130, 618];
m = length(integrationTimes);
measured = cell(m,1);
xPoints = [100, 500, 900, 1200];
yPoints = [100, 500, 900];
for i = 1:m
    setSetting('integrationTime', integrationTimes(i));
    setSetting('saveFolder', fullfile(experiment, 'whiteComparison'));
    readHSIData('whiteReflectance', target, experiment);
    fileConditions = getFileConditions('whiteReflectance', target);
    [measured{i}, curveNames] = getRepresentativePoints(fileConditions, xPoints, yPoints);
end 

for i = 1:m
    measured{i} = reshape(measured{i}, [3*4 401]);
end 
expected = measured{m};
nmses = zeros(13, m-1);
for i = 1:(m-1)
    nmses(1:12,i) = applyRowFunc(@nmse, measured{i}, expected);
    nmses(13,i) = mean(nmses(1:12,i));
end 
nmsesTable = array2table(nmses,'VariableNames',{'NMSE(618vs1380)','NMSE(618vs800)','NMSE(618vs2130)'});
nmsesTable.Names = [curveNames; 'Average'];

rmses = zeros(13, m-1);
for i = 1:(m-1)
    rmses(1:12,i) = applyRowFunc(@rmse, measured{i}, expected);
    rmses(13,i) = mean(rmses(1:12,i));
end 
rmsesTable = array2table(rmses,'VariableNames',{'RMSE(618vs1380)','RMSE(618vs800)','RMSE(618vs2130)'});
rmsesTable.Names = [curveNames; 'Average'];
save( fullfile(getSetting('savedir'), getSetting('saveFolder'), 'errorTables.mat'), 'rmsesTable', 'nmsesTable');

%% Colorchart middle with different normalizations
experiment = strcat('testCalibration', version);
setSetting('experiment', experiment);
setSetting('isRotated', false);
setSetting('integrationTime', 618);
target = strcat('colorchart', 'Middle');
readHSIData('colorchart', target, experiment);
allowRoiSelection = true;
setSetting('colorPatchOrder', 'redRight');

setSetting('saveFolder', fullfile(experiment, 'raw'));
setSetting('normalization', 'raw');
evaluateColorchart(target, allowRoiSelection); 
    
normalizations = { 'bandmax', 'uniSpectrum', 'byPixel'};
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
    readHSIData('colorchart', target, experiment);
    [tables{k}, measuredSpectra{k}, adjustedSpectra{k}, alphas{k}] = evaluateColorchart(target, allowRoiSelection); 
end 

[tables] = exportSimilarityTables(tables, measuredSpectra, adjustedSpectra, alphas, confs);

endRun;
function [] = importImagesToMat(dataDate, experiment, targetInfo, configuration, database)
%% IMPORTIMAGESTOMAT prepares data by reading hsm files from the db as .mat files 

if nargin < 4
    configuration = 'singleLightClose';
end 
if nargin < 5 
    dataBase = 'calib';
end

normalization = 'byPixel';
initialization;

targets = targetInfo.Targets; 
integrationTimes = targetInfo.integrationTimes;
types = targetInfo.Types;

for i = 1:numel(targetInfo)       
    setSetting('integrationTime', integrationTimes(i));
    readHSIData(types{i}, targets{i}, experiment);
    setSetting('isRotated', false);
end 

end 
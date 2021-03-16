function [searchName, integrationTime] = getFilename(configuration, content, integrationTime)

%% getFilename Gets the respective filename for configuration value
%   'configuration' and content value 'content' for the case of hsi
%   calibration
warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
dataTable = readtable(fullfile(getSetting('datasetSettingsDir'), 'dataSetCharacteristics.csv'));
if nargin < 3
    setId = strcmp(dataTable.configuration, configuration) & strcmp(dataTable.content, content);
else
    setId = strcmp(dataTable.configuration, configuration) & strcmp(dataTable.content, content) & (dataTable.exposureTime == integrationTime);
end

if sum(setId)  > 1
    dataDate = getSetting('dataDate');
    setId = (dataTable.date == str2num(dataDate)) & setId; 
end

searchName = dataTable.filename{setId};

if nargin > 2
    if integrationTime ~= dataTable.exposureTime(setId)
        warning('Integration time in the settings and in the retrieved file differs.');
    end
end

integrationTime = dataTable.exposureTime(setId);
setSetting('integrationTime', integrationTime);

end
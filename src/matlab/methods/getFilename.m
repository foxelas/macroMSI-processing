function searchName = getFilename(configuration, content)
%% getFilename Gets the respective filename for configuration value
%   'configuration' and content value 'content' for the case of hsi 
%   calibration
    warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
    dataTable = readtable( fullfile(getSetting('datasetSettingsDir'), 'dataSetCharacteristics.csv')); 
    setId = strcmp(dataTable.configuration, configuration) & strcmp(dataTable.content, content);
    searchName = dataTable.filename{setId};
end 
function searchVal = getIndexValue(searchName, dataDate)
%% getIndexValue Gets the respective index for value name 'searchName'
%   and case name 'dataDate' for the case of hsi calibration

    datasetCharacteristics = delimread( fullfile(getSetting('datasetSettingsDir'), 'dataSetCharacteristics.csv'), ',', {'mixed'}).mixed;
    setId = find(strcmp(datasetCharacteristics(2:end,1), searchName)  &  ([datasetCharacteristics{2:end,2}] == str2num(dataDate))');
    searchVal = datasetCharacteristics(setId + 1, 3);
    searchVal = searchVal{1};
end 
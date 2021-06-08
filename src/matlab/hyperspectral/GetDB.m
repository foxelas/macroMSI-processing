function dataTable = GetDB()
%%Returns the db structure
dataTable = readtable(fullfile(getSetting('datasetSettingsDir'), strcat(getSetting('database'), 'DB.xlsx')), 'Sheet','capturedData');
end 

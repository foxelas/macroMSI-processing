function [filename, tableId, outRow] = getFilename(configuration, content, integrationTime, target, dataDate, id)
%% getFilename Gets the respective filename for configuration value
%   arguments are received in the order of 
%     'configuration' [light source]
%     'content' [type of catpured object]
%     'integrationTime' [value of integration time]
%     'target' [details about captured object]
%     'dataDate' [catpureDate]
%     'id' [number value for id ]

warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
dataTable = readtable(fullfile(getSetting('datasetSettingsDir'), strcat(getSetting('database'), 'DB.xlsx')), 'Sheet','capturedData');

if ~isempty(integrationTime)
    initialIntegrationTime = integrationTime;
else 
    initialIntegrationTime = [];
end 

setId = ismember(dataTable.Configuration, configuration);
if nargin >= 2 && ~isempty(content)
    setId = setId & ismember(dataTable.Content, content);
end 
if nargin >= 3 && ~isempty(integrationTime)
    setId = setId & ismember(dataTable.IntegrationTime, integrationTime);
end 
if nargin >= 4 && ~isempty(target)
    setId = setId & ismember(dataTable.Target, target);
end 
if nargin >= 5 && ~isempty(dataDate)
    setId = setId & ismember(dataTable.CaptureDate, str2num(dataDate));
end 
if nargin >= 6 && ~isempty(id)
    setId = setId & ismember(dataTable.ID, id);
end 

outRow = dataTable(setId,:);
if sum(setId) > 1 
    warning('Taking the first from multiple rows that satisfy the conditions.');
    outRow = outRow(1,:);
end 

filename = outRow.Filename{1};
tableId = outRow.ID;

if outRow.IntegrationTime ~= initialIntegrationTime
    setSetting('integrationTime', integrationTime);
end 

if nargin >= 3 && ~isempty(integrationTime) &&  integrationTime ~= outRow.IntegrationTime
    warning('Integration time in the settings and in the retrieved file differs.');
%     setSetting('integrationTime', integrationTime);
end


end
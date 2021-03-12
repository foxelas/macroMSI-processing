%% Setup 
experiment = 'testStocmach';
dataDate = '20201218';
%configuration = 'singleLightClose';
integrationTime = 200;
initialization;

configurations = { 'doubleLightClose'};
n = numel(configurations);

for i = 1:n
    configuration = configurations{i};
    setSetting('configuration', configuration)
    setSetting('saveFolder', fullfile(experiment, configuration));
    readWhite(dataDate, integrationTime, true, false, experiment, configuration, []);
end
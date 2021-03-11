%% Setup 
experiment = 'forEMBC';
dataDate = '20201218';
%configuration = 'singleLightClose';
integrationTime = 200;
initialization;

configurations = { 'singleLightClose', 'doubleLightClose'};
n = numel(configurations);

%% Read white
readBlackWhite = false;
if readBlackWhite
    for i = 1:n
        configuration = configurations{i};
        setSetting('configuration', configuration)
        setSetting('saveFolder', fullfile(experiment, configuration));
        readWhite(dataDate, integrationTime, true, false, experiment, configuration, []);
    end
end


%% Compare light setting 
tables = cell(n,1);
measuredSpectra = cell(n,1);
adjustedSpectra = cell(n,1);
allowRoiSelection = true; 

for i = 1:n
    configuration = configurations{i};
    setSetting('configuration', configuration)
    setSetting('saveFolder', fullfile(experiment, configuration));
    %plotSaveDir = fullfile(getSetting('savedir'), getSetting('saveFolder'));

    getRepresentativePoints('whiteReflectance');

    [tables{i}, measuredSpectra{i}, adjustedSpectra{i}] = evaluateColorchart('colorchart', allowRoiSelection);
 
end
%% Compare filter setting 
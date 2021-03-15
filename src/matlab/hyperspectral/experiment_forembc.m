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

confs = cell(n,2);
normalization = getSetting('normalization');
for i = 1:n
    configuration = configurations{i};
    setSetting('configuration', configuration)
    setSetting('saveFolder', fullfile(experiment, configuration));

%     getRepresentativePoints('whiteReflectance');
    confs(i, 1:2) = deal({configuration, normalization});
    [tables{i}, measuredSpectra{i}, adjustedSpectra{i}] = evaluateColorchart('colorchart', allowRoiSelection); 
end

 
[tables] = exportSimilarityTables(tables, measuredSpectra, adjustedSpectra, confs);


%% Compare filter setting 
%% Setup 
startRun;
experiment = 'forEMBC';
dataDate = '20201218';
%configuration = 'singleLightClose';
integrationTime = 200;
normalization = 'byPixel';

initialization;

configurations = { 'singleLightClose', 'doubleLightClose'};
normalizations = {'bandmax', 'uniSpectrum', 'byPixel'};
n = numel(configurations);
v = numel(normalizations);

%% Read white
readBlackWhite = false;
if readBlackWhite
    for i = 1:n
        configuration = configurations{i};
        setSetting('configuration', configuration)
        setSetting('saveFolder', fullfile(experiment, configuration));
        readWhite(dataDate, integrationTime, experiment, configuration, []);
    end
end


%% Compare light setting
m = n*v;
tables = cell(m,1);
measuredSpectra = cell(m,1);
adjustedSpectra = cell(m,1);
alphas = cell(m,1);
allowRoiSelection = false; 

confs = cell(m,2);
normalization = getSetting('normalization');
k = 0; 
for i = 1:n
    for j = 1:v
        k = k + 1;
        configuration = configurations{i};
        setSetting('configuration', configuration)
        normalization = normalizations{j};
        setSetting('normalization', normalization);
        setSetting('saveFolder', fullfile(experiment, strcat(configuration, '_', normalization)));

    %     getRepresentativePoints('whiteReflectance');
        confs(k, 1:2) = deal({configuration, normalization});
        [tables{k}, measuredSpectra{k}, adjustedSpectra{k}, alphas{k}] = evaluateColorchart('colorchart', allowRoiSelection); 
    end 
end

[tables] = exportSimilarityTables(tables, measuredSpectra, adjustedSpectra, alphas, confs);

endRun;

%% Compare filter setting 
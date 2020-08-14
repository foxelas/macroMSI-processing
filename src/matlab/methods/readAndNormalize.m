function [preprocMsi, whiteReference, specimenMask] = readAndNormalize(i, ID, msiType, removeBg, normType)

%% Load Data
groupId = i;
poiId = groupId2poiId(groupId, ID); %id = ID(78);%fixed data
infile = fullfile(getSetting('systemdir'), 'infiles', strcat('group_', num2str(groupId), '.mat'));
load(infile, 'raw', 'specimenMask');

close all;

%% Print MSI
[msi, whiteReference, specimenMask, ~, ~, ~] = getImage(groupId, msiType, removeBg, false, 'none');
setSetting('saveImages', true);
setSetting('plotName', fullfile(getSetting('savedir'), getSetting('common'), strcat(num2str(groupId), '_', 'msi.png')));
plotFunWrapper(1, @plotMSI, msi);

%% Normalize MSI
[msiNorm] = getImage(groupId, msiType, removeBg, false, normType);

setSetting('saveImages', true);
setSetting('plotName', fullfile(getSetting('savedir'), getSetting('common'), strcat(num2str(groupId), '_', 'msiNorm.png')));
plotFunWrapper(2, @plotMSI, msiNorm);

preprocMsi = msiNorm;

end

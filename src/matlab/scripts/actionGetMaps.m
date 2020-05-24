% k = 51;
% infile = fullfile(getSetting('systemdir'), 'infiles', strcat('poi_', num2str(k), '.mat'));
% load(infile, 'poiName', 'poiRAW', 'poiSegmentMask', ...
%         'roiSeeds', 'measuredSpectrum', 'poiWhite');
k = 14; %fixed data 
k = 13; %unfixed data
infile = fullfile(getSetting('systemdir'), 'infiles', strcat('group_', num2str(k), '.mat'));
load(infile, 'raw', 'specimenMask');

[map] = getMap(raw, 'opticalDensityMelanin', specimenMask);
[map] = getMap(raw, 'opticalDensityHemoglobin', specimenMask);
[map] = getMap(raw, 'hemoglobin', specimenMask);
[map] = getMap(raw, 'totalMelanin', specimenMask);

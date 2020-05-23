k = 51;
infile = fullfile(getSetting('systemdir'), 'infiles', strcat('poi_', num2str(k), '.mat'));
load(infile, 'poiName', 'poiRAW', 'poiSegmentMask', ...
        'roiSeeds', 'measuredSpectrum', 'poiWhite');
    
[map] = getMap(poiRAW, 'totalMelanin');


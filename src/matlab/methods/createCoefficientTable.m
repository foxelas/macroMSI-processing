function [] = createCoefficientTable(ID, CompleteSpectra)
%     CREATECOEFFICIENTTABLE creats the coefficient table from
%     the measured spectra.
%
%     Usage:
%     createCoefficientTable(ID, CompleteSpectra)

disp('Reading coefficients from excel file...')
pixelValueSelectionMethods = {'green', 'rms', 'adjusted'};
Coefficients = ones(length(ID), 3, 7);
xls = dir(fullfile(getSetting('systemdir'), 'coeff*'));
src = fullfile(getSetting('systemdir'), xls(1).name);
msiN = length(ID);
for k = 1:msiN
    idd = ID(k);
    load(fullfile(getSetting('systemdir'), 'infiles', strcat('poi_', num2str(idd.Index), '.mat')), ...
        'poiRAW', 'roiSeeds');
    % Retrieve MSI data
    raw = poiRAW;
    % Retrieve spectrum data
    refmeasured = CompleteSpectra(k, :);

    for m = 1:length(pixelValueSelectionMethods)
        msi = raw2msi(raw, pixelValueSelectionMethods{m});
        x = roiSeeds(1);
        y = roiSeeds(2);
        refmsi = im2uint16(squeeze(msi(:, y, x)))';
        xlswrite2(src, refmeasured, 4, 'B3');
        xlswrite2(src, refmsi, 4, 'D407:J407');
        ctemp = xlsread(src, 4, 'D414:J414');
        Coefficients(k, m, :) = ctemp;
    end
end
save(fullfile(getSetting('systemdir'), 'precomputedParams.mat'), 'Coefficients', '-append');
end

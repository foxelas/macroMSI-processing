close all; clc;
k = 20;
msiType = 'extended'; %'extended'; % 'max';

% for all pixels without bg removal
%[msi, ~, ~, ~, ~, ~] = getImage(k, msiType, true);
%figure(4);imshow(squeeze(msi(1,:,:)));
%[msi] = getNormMSI(msi);
m = 150;
n = 150;
c = size(msi, 1);

patchNames = {'Healthy patch', 'Atypical patch', 'Cancer patch'};
patches = {msi(:, 1187:(1187 + n - 1), 438:(438 + m - 1)), msi(:, 820:(820 + n - 1), 685:(685 + m - 1)), msi(:, 699:(699 + n - 1), 444:(444 + m - 1))};

setSetting('saveImage', false);
nbins = 10;
for i = 1:length(patchNames)
    patch = patches{i};
    c = getBandCorrelation(patch);
    plotBandCorrelations(c, msiType, sprintf('Band Correlation for %s', patchNames{i}), i);
end
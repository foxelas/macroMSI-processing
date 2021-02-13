close all; clc;
k = 20;
msiType = 'extended'; %'extended'; % 'max';

% for all pixels without bg removal
removebg = false;
[msi, ~, ~, ~, ~, ~] = getImage(k, msiType, true);
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
    plots(i, @plotMSI, patch);
    title(patchNames{i});
    spectralHist = zeros(c, nbins);
    for j = 1:c
        spectralHist(j, :) = histcounts(squeeze(patch(j, :, :)), nbins);
    end
    spectralHists{i} = spectralHist;
end

simHealthyAtypical = getSpectralHistogramSimilarity(spectralHists{1}, spectralHists{2})
simHealthyCancer = getSpectralHistogramSimilarity(spectralHists{1}, spectralHists{3})
simCancerAtypical = getSpectralHistogramSimilarity(spectralHists{3}, spectralHists{2})

simHealthyHealthy = getSpectralHistogramSimilarity(spectralHists{1}, spectralHists{1})
simAtypicalAtypical = getSpectralHistogramSimilarity(spectralHists{2}, spectralHists{2})
simCancerCancer = getSpectralHistogramSimilarity(spectralHists{3}, spectralHists{3})

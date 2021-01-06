function [reorderedSpectra, labels] = reorderSpectra(target, chartColorOrder, spectraColorOrder, wavelengths, spectralWavelengths)
% match chartColorOrder according to spectralColorOrder 
% i.e. match babel order to colorchart order 
[~, idx] = ismember(spectraColorOrder, chartColorOrder);
idx = nonzeros(idx); 
[~, idx2] = ismember(spectralWavelengths', wavelengths);

targetDecim = target(:, idx2);
reorderedSpectra = targetDecim(idx, :);

labels = spectraColorOrder;
end
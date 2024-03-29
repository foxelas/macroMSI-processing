function [reorderedSpectra, labels] = reorderSpectra(target, chartColorOrder, spectraColorOrder, wavelengths, spectralWavelengths)
% match chartColorOrder according to spectralColorOrder
% i.e. match babel order to colorchart order

if size(target, 2) == 401
    wavelengths = wavelengths;
elseif size(target, 2) == 161
    wavelengths = [380:540]';
else
    wavelengths = [541:780]';
end

[~, idx] = ismember(spectraColorOrder, chartColorOrder);
idx = nonzeros(idx);
[~, idx2] = ismember(spectralWavelengths', wavelengths);

idx2 = nonzeros(idx2);

targetDecim = target(:, idx2);
reorderedSpectra = targetDecim(idx, :);

labels = spectraColorOrder;
end
close all;

searchdir = getSetting('datadir');
datalist = dir(fullfile(searchdir, '**', 'white.csv'));
n = length(datalist);
spectrum = zeros(n, 401);
figure(1)
hold on
for i = 1:n
    spectrum(i, :) = readSpectrum(fullfile(datalist(i).folder, datalist(i).name))';
    plot(380:780, spectrum(i, :))
end
referenceSpectrum = mean(spectrum);
referenceSpectrum401 = referenceSpectrum;
referenceSpectrum81 = decimate(referenceSpectrum, 5);
h = plot(380:780, referenceSpectrum, 'LineWidth', 5, 'DisplayName', 'Average');
legend(h);
hold off
xlabel('wavelength (nm)');
ylabel('reflectance (a.u.)')
title('White Reference Spectra')

save(fullfile('parameters', 'referenceSpectrum.mat'), 'referenceSpectrum401', 'referenceSpectrum81');
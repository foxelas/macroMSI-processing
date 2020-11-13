function[] = plotSpectra(spectra, wavelengths, names, figTitle, fig)

if isempty(wavelengths)
    wavelengths = 380:5;780; 
end 

if isempty(names)
    names = [];
end

if ~iscell(names)
    names = {names};
end 

if isempty(figTitle)
    figTitle = 'Calculated Spectra'; 
end 

lineColorMap = getLineColorMap('custom', names);
key = keys(lineColorMap);

hold on
for i = 1:length(names)
    h = plot(wavelengths, spectra(i,:), 'DisplayName', key{i}, 'Color', lineColorMap(key{i}));
end
hold off

legend(h, 'Location', 'northwest', 'FontSize', 15)
xlabel('Wavelength (nm)', 'FontSize', 15);
ylabel('Reflectance (a.u.)', 'FontSize', 15);
title (figTitle)

savePlot(fig);

end
%% Configuration plots about illumination, sensitivity etc
options.saveOptions.plotName = strcat(options.saveOptions.savedir, 'illuminationPlot.jpg');
plots('illumination', 1, [], '', 'wavelength', wavelength, 'illumination', illumination, 'saveOptions', options.saveOptions);
options.saveOptions.plotName = strcat(options.saveOptions.savedir, 'sensitivityPlot.jpg');
plots('sensitivity', 2, [], '', 'wavelength', wavelength, 'sensitivity', sensitivity, 'saveOptions', options.saveOptions);
options.saveOptions.plotName = strcat(options.saveOptions.savedir, 'illuminationAndSensitivityPlot.jpg');
plots('illuminationAndSensitivity', 3, [], '', 'wavelength', wavelength, 'illumination', illumination, 'sensitivity', sensitivity, 'saveOptions', options.saveOptions);
% end of Configuration plots about illumination, sensitivity etc
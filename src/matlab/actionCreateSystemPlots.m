%% Configuration plots about illumination, sensitivity etc
%options.saveOptions.BW = false;
parts = regexp(options.saveOptions.savedir, '\', 'split'); 
parts(end) =  [];
saveDir = fullfile(strjoin(parts, '\'), 'general');
options.saveOptions.plotName = fullfile(saveDir, 'illuminationPlot');
plots('illumination', 1, [], '', 'Wavelength', wavelength, 'Illumination', illumination, 'SaveOptions', options.saveOptions);
options.saveOptions.plotName = fullfile(saveDir, 'sensitivityPlot'); 
plots('normSensitivity', 2, [], '', 'Wavelength', wavelength, 'Sensitivity', sensitivity, 'SaveOptions', options.saveOptions);
options.saveOptions.plotName = fullfile(saveDir, 'illuminationAndSensitivityPlot');
plots('illuminationAndSensitivity', 3, [], '', 'Wavelength', wavelength, 'Illumination', illumination, 'Sensitivity', sensitivity, 'SaveOptions', options.saveOptions);
% end of Configuration plots about illumination, sensitivity etc

%% Count and analyze the contents of the dataset
datasetBreakdown(ID, options);
%         Breakdown of the dataset:
%         Unfixed: 43, Fixed: 45, Cut: 39
%         Normal: 74, Cancerous: 53
%
%         Among unfixed data:
%         Normal: 23, Cancerous: 20

% end of  Count and analyze the contents of the dataset
%% Configuration plots about illumination, sensitivity etc
options.saveOptions.plotName = strcat(options.saveOptions.savedir, 'illuminationPlot.jpg');
plots('illumination', 1, [], '', 'Wavelength', wavelength, 'Illumination', illumination, 'SaveOptions', options.saveOptions);
options.saveOptions.plotName = strcat(options.saveOptions.savedir, 'sensitivityPlot.jpg');
plots('sensitivity', 2, [], '', 'Wavelength', wavelength, 'Sensitivity', sensitivity, 'SaveOptions', options.saveOptions);
options.saveOptions.plotName = strcat(options.saveOptions.savedir, 'illuminationAndSensitivityPlot.jpg');
plots('illuminationAndSensitivity', 3, [], '', 'Wavelength', wavelength, 'Illumination', illumination, 'Sensitivity', sensitivity, 'SaveOptions', options.saveOptions);
% end of Configuration plots about illumination, sensitivity etc

%% Count and analyze the contents of the dataset
[~, idx, ~] = unique(strcat({ID.Csvid}, {ID.T}));
measuredSpectraCount = length(idx);
normalFixedCount = 0;
normalCutCount = 0;
normalUnfixedCount = 0;
fixedCount = 0;
cutCount = 0;
unfixedCount = 0;
for i = 1:length(idx)
    if ~(ID(idx(i)).IsFixed)
        unfixedCount = unfixedCount + 1;
        if ID(idx(i)).IsNormal
            normalUnfixedCount = normalUnfixedCount + 1;
        end
    elseif (ID(idx(i)).IsFixed && ~(ID(idx(i)).IsCut))
        fixedCount = fixedCount + 1;
        if ID(idx(i)).IsNormal
            normalFixedCount = normalFixedCount + 1;
        end
    else
        cutCount = cutCount + 1;
        if ID(idx(i)).IsNormal
            normalCutCount = normalCutCount + 1;
        end
    end
end
normalCount = normalUnfixedCount + normalFixedCount + normalCutCount;
cancerCount = length(idx) - normalCount;
cancerUnfixedCount = unfixedCount - normalUnfixedCount;

dataCount = strcat('Breakdown of the dataset:\nTotal: ', num2str(length(idx)), ...
    ', Unfixed: ', num2str(unfixedCount), ', Fixed: ', num2str(fixedCount), ...
    ', Cut: ', num2str(cutCount), '\n', 'Normal: ', num2str(normalCount), ...
    ', Cancerous: ', num2str(cancerCount), '\n\nAmong unfixed data:\n', ...
    'Normal: ', num2str(normalUnfixedCount), ', Cancerous: ', num2str(unfixedCount-normalUnfixedCount),'\n',...
    'Among fixed data:\nNormal: ', num2str(normalFixedCount), ', Cancerous: ', num2str(fixedCount-normalFixedCount), '\n\n', ...
    'Among cut data:\nNormal: ', num2str(normalCutCount), ', Cancerous: ', num2str(cutCount-normalFixedCount), '\n\n');

disp(dataCount);

fileID = fopen( fullfile( options.savedir, 'data_count.txt'), 'w');
fprintf(fileID, outputLog);
fclose(fileID);

%         Breakdown of the dataset:
%         Unfixed: 43, Fixed: 45, Cut: 39
%         Normal: 74, Cancerous: 53
%
%         Among unfixed data:
%         Normal: 23, Cancerous: 20

% end of  Count and analyze the contents of the dataset
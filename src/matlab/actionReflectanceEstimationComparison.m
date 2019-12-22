%Start of reflectance estimation comparison
options.action = 'Refest';

w = warning('off', 'all');

methods = {'MSI-Simple', 'MSI-SpatioSpectral', 'MSI-Spatial', 'RGB-Simple'};
methodsN = length(methods);
multipleReconstructions = zeros(msiN, length(wavelength), methodsN);
nmses = zeros(msiN, methodsN);
gfcs = zeros(msiN, methodsN);

for k = 1:msiN
    % Retrieve MSI data

    infile = fullfile(options.systemdir, 'infiles', strcat('poi_', num2str(k), '.mat'));
    load(infile, 'poiName', 'poiRAW', 'poiSegmentMask', ...
        'roiSeeds', 'measuredSpectrum', 'poiWhite', 'sigma2');

    msi = poiRAW;
    mask = poiSegmentMask;
    measured = measuredSpectrum;
    rgb = poiWhite;
    %options.sigma2 = sigma2;

    [reconstructedArray, gfcArray, nmseArray] = getMultipleReflectanceReconstructions(msi, rgb, mask, measured, ID(k), options, methods);
    multipleReconstructions(k, :, :) = reconstructedArray;
    nmses(k, :) = nmseArray;
    gfcs(k, :) = gfcArray;

end
warning(w);

%% Export results
nmseTable = array2table(nmses, ...
    'VariableNames', cellfun(@(x) strrep(x, '-', ''), methods, 'un', 0));
gfcTable = array2table(gfcs, ...
    'VariableNames', cellfun(@(x) strrep(x, '-', ''), methods, 'un', 0));
errorInfo = GetErrorInfoStruct(nmses, gfcs);
saveOptions = options.saveOptions;
saveOptions.saveImages = true;
saveOptions.showImages = true;
for i = 1:methodsN
    nmseBars = getNmseForBarPlot(nmses(:, i), ID, methods{i}, saveOptions);
    saveOptions.plotName = fullfile(saveOptions.savedir, outputFolderMap('reflectanceEstimationPerformance'), strcat('hist_', methods{i}));
    plotGFCHistogram(gfcs(:, i), 1, saveOptions);
end

msiId = find(strcmp(methods, 'MSI-Simple'));
EstimatedSpectra = multipleReconstructions(:, :, msiId);
rgbId = find(strcmp(methods, 'RGB-Simple'));
EstimatedRGBSpectra = multipleReconstructions(:, :, rgbId);
filename = mkdir_custom(fullfile(saveOptions.savedir, outputFolderMap('features'), 'out.mat'));
if exist(filename, 'file')
    save(filename, 'EstimatedSpectra', 'EstimatedRGBSpectra', 'nmseTable', 'gfcTable', 'errorInfo', 'multipleReconstructions', 'methods', '-append');
else
    save(filename, 'EstimatedSpectra', 'EstimatedRGBSpectra', 'nmseTable', 'gfcTable', 'errorInfo', 'multipleReconstructions', 'methods');
end
%End of reflectance estimation comparison


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function errorInfo = GetErrorInfoStruct(nmse, gfc)
nmse = nmse';
gfc = gfc';
errorInfo = struct('avgnmse', num2cell(mean(nmse, 2)), 'minnmse', num2cell(min(nmse, [], 2)), ...
    'maxnmse', num2cell(max(nmse, [], 2)), 'stdnmse', num2cell(std(nmse, [], 2)), 'stdeSEM', ...
    num2cell(std(nmse, [], 2)./size(nmse, 2)), 'avggfc', num2cell(mean(gfc, 2)), ...
    'mingfc', num2cell(min(gfc, [], 2)), 'maxgfc', num2cell(max(gfc, [], 2)), 'stdgfc', ...
    num2cell(std(gfc, [], 2)));
errorInfo = struct2table(errorInfo);
end


function nmseBars = getNmseForBarPlot(nmses, ID, method, saveOptions)
i = 0;
nmseBars = zeros(6, 1);
for state = {'unfixed', 'fixed', 'cut'}
    for malignancy = 0:1
        i = i + 1;
        idx = [ID.IsBenign] == malignancy & strcmp({ID.Type}, state);
        nmseBars(i) = sum(nmses(idx)) / length(ID);
    end
    %fprintf('NRMSE for %s = %.4f\n', state{1}, nrmse);
end
nmseBars = reshape(nmseBars, [2, 3]);

saveOptions.plotName = fullfile(saveOptions.savedir, outputFolderMap('reflectanceEstimationPerformance'), method);
plotReconstructionPerformanceBars(nmseBars, {'malignant', 'benign'}, '', 1, saveOptions);

fprintf('NRMSE overall = %.4f\n', sum(nmseBars(:)));

end

function [] = plotVisualResultAfterClassification(cancerProb, predictedLabels, coordinates, orderedRois, segmentMasks, sRGB, outName, options, fig, imgTitle)
    
    key = {0, 1, 2}; 
    value = {'Benign', 'Malignant', 'Atypical'};
    labelMap = containers.Map(key, value);
    predictedClasses = values(labelMap, num2cell(predictedLabels));
    
    
    [M, N, ~] = size(sRGB);
    malignancyMap = zeros(M,N);
    for i = orderedRois'
        malignancyMap(logical(segmentMasks{i})) = cancerProb(i);
    end
    
    options.saveOptions.plotName = fullfile(options.saveOptions.savedir, outputFolderMap('visualTool'), outName);
    %labels = {isPositive, ~[ID(roiIndexes).IsBenign]};
    plotVisualResult(sRGB, malignancyMap, imgTitle, predictedClasses, coordinates, 'hsv', false, fig, options.saveOptions)
end
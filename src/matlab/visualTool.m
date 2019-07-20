for k = [20, 19, 18]
    infile = fullfile(options.systemdir, 'infiles', strcat('group_', num2str(k), '.mat'));
    load(infile, 'raw', 'whiteReference', 'specimenMask');
    figure(5); imshow(whiteReference);
    z = find([ID.Group] == k, 1);
    sRGB = createSRGB(raw, 'medium', ID(z), options, 'cmccat2000', specimenMask);
    

    trueLabels = [1 1 1 1 0 0];
    options.saveOptions.saveImages = true;
    roiIndexes = find([ID.Group] == k);
    rois = length(roiIndexes);
    coordinates = [[ID(roiIndexes).Originx]; [ID(roiIndexes).Originy]]';
    segmentMasks = cell(length(roiIndexes), 1);
    for i = 1:length(roiIndexes) 
        infile = fullfile(options.systemdir, 'infiles', strcat('poi_', num2str(roiIndexes(i)), '.mat'));
        load(infile, 'segmentMask');
        segmentMasks{i} = segmentMask;
    end
    [~, orderedRois] = sort(cellfun(@(x) sum(x(:)), segmentMasks), 'descend');
    
    if k == 20
        cancerProb = 1 - [0.40334867 0.41269009 0.33632419 0.31629984 0.29190125 0.35072469];
        predictedLabels = [1 1 0 0 1 0];
        getImageForClassifier(cancerProb, predictedLabels, coordinates, orderedRois, segmentMasks, sRGB, 'Knn_unfixed', options, 1, 'Unfixed');

        cancerProb = 1 - [0.615 0.605 0.56  0.49  0.655 0.53];
        predictedLabels = [1 1 1 0 1 1];
        getImageForClassifier(cancerProb, predictedLabels, coordinates, orderedRois, segmentMasks, sRGB, 'RF_unfixed', options, 2, 'Unfixed');
    elseif k == 19
        cancerProb = 1 - [0.33009663 0.35727159 0.42860166 0.40599089  0.38940433 0.34980716 0.35152617 0.3277389];
        predictedLabels = [1 1 1 1 1 1 1 1];
        getImageForClassifier(cancerProb, predictedLabels, coordinates, orderedRois, segmentMasks, sRGB, 'Knn_fixed', options, 1, 'Fixed');

        cancerProb = 1 - [0.46  0.545 0.545 0.6  0.525 0.555 0.605 0.61 ];
        predictedLabels = [ 0 1 1 1 1 1 1 1];
        getImageForClassifier(cancerProb, predictedLabels, coordinates, orderedRois, segmentMasks, sRGB, 'RF_fixed', options, 2, 'Fixed');
    elseif k == 18
        cancerProb = 1 - [0.34954588 0.38105323 0.33074516 0.31270676 0.30913618 0.3470074  0.35892913 0.32471956];
        predictedLabels = [1 1 1 1 1 1 0 0];
        getImageForClassifier(cancerProb, predictedLabels, coordinates, orderedRois, segmentMasks, sRGB, 'Knn_sectioned', options, 1, 'Sectioned');

        cancerProb = 1 - [0.51  0.585 0.61  0.615 0.59  0.63  0.49  0.475];
        predictedLabels = [1 1 1 1 1 1 0 0];
        getImageForClassifier(cancerProb, predictedLabels, coordinates, orderedRois, segmentMasks, sRGB, 'RF_sectioned', options, 2, 'Sectioned');
    end
     
end

function [] = getImageForClassifier(cancerProb, predictedLabels, coordinates, orderedRois, segmentMasks, sRGB, outName, options, fig, imgTitle)
    [M, N, ~] = size(sRGB);
    malignancyMap = zeros(M,N);
    for i = orderedRois'
        malignancyMap(logical(segmentMasks{i})) = cancerProb(i);
    end
    
    options.saveOptions.plotName = fullfile(options.saveOptions.savedir, '11-Visual Tool', outName);
    %labels = {isPositive, ~[ID(roiIndexes).IsBenign]};
    plotVisualResult(sRGB, malignancyMap, imgTitle, predictedLabels, coordinates, 'hsv', false, fig, options.saveOptions)
end
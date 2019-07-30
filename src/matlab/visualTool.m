for k = [20, 19, 18]

    outfile = fullfile(options.saveOptions.savedir, '10-sRGB', strcat('group_', num2str(k), '.mat'));
    load(outfile, 'sRGB');

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
%         %SVM-rbf-3.0-auto-True-harsh_penalty-spect+mlbp-ICA-20-ICA-20
%         % 3scales
%         cancerProb =  [ 0.79466855 0.63459283  0.76740779 0.73678306 0.52148642 0.72221257 ];
%         predictedLabels = [1 1 1 1 0 1];
        cancerProb = [0.83883908 0.84836129  0.70056079 0.67005854 0.79753921 0.70401719];
        predictedLabels = [1 1 1 1 1 1];
        getImageForClassifier(cancerProb, predictedLabels, coordinates, orderedRois, segmentMasks, sRGB, 'Svm_unfixed', options, 3, 'Unfixed');
        
        % KNN-3-chebyshev-distance-spect+mlbp-PCA-20-None-None
        % 2 scales
        cancerProb =  [1.         0.6845067  0.64820811 0.3260943  0.31588529 0.       ];
        predictedLabels = [1 1 1 0 0 0];
        getImageForClassifier(cancerProb, predictedLabels, coordinates, orderedRois, segmentMasks, sRGB, 'Knn_unfixed', options, 1, 'Unfixed');
        %Random Forest-100-gini-log2-none-balanced-spect+mlbp-ICA-20-ICA-20
        %2 scales
        cancerProb =  [0.55 0.56 0.56 0.54 0.61 0.45];
        predictedLabels = [1 1 1 1 1 0];
        getImageForClassifier(cancerProb, predictedLabels, coordinates, orderedRois, segmentMasks, sRGB, 'RF_unfixed', options, 2, 'Unfixed');
        
        trueLabels = [1 1 1 1 0 0];
        getImageForClassifier([0 0 0 0 0 0], trueLabels, coordinates, orderedRois, segmentMasks, sRGB, 'GD_unfixed', options, 2, 'Fixed');

    elseif k == 19
%         cancerProb =  [0.6683761  0.78546708 0.85105244 0.83916144  0.77098218 0.76302766 0.80217468 0.78879581 ];
%         predictedLabels = [ 1 1 1 1 1 1 1 1 ];
        cancerProb = [0.74922035 0.79203226 0.77949437 0.81355963  0.74810811 0.7854452  0.80657934 0.86976086 ];
        predictedLabels = [ 1 1 1 1 1 1 1 1 ];
        getImageForClassifier(cancerProb, predictedLabels, coordinates, orderedRois, segmentMasks, sRGB, 'Svm_fixed', options, 3, 'Fixed');
        
        cancerProb =  [1.         1.         1.         1.   1.         1.         0.32505164 0.32528109];
        predictedLabels = [1 1 1 1 1 1 0 0 ];
        getImageForClassifier(cancerProb, predictedLabels, coordinates, orderedRois, segmentMasks, sRGB, 'Knn_fixed', options, 1, 'Fixed');

        cancerProb = [ 0.56 0.62 0.57 0.64 0.58 0.54 0.58 0.64];
        predictedLabels = [ 1 1 1 1 1 1 1 1 ];
        getImageForClassifier(cancerProb, predictedLabels, coordinates, orderedRois, segmentMasks, sRGB, 'RF_fixed', options, 2, 'Fixed');
        
        trueLabels = [1 1 1 1 1 1 0 0 ];
        getImageForClassifier([0 0 0 0 0 0 0 0], trueLabels, coordinates, orderedRois, segmentMasks, sRGB, 'GD_fixed', options, 2, 'Fixed');

    elseif k == 18
%         cancerProb =  [0.82120848 0.67357772 0.85494806 0.62093993 0.711827   0.76968064  0.67605464 0.54017125  ];
%         predictedLabels = [1 1 1 1 1 1 1 0 ];
        cancerProb = [0.69209988 0.52010556 0.77349548 0.60031225 0.68874827 0.72247141  0.3476265  0.30541984  ];
        predictedLabels = [1 0 1 1 1 1 0 0  ];
        getImageForClassifier(cancerProb, predictedLabels, coordinates, orderedRois, segmentMasks, sRGB, 'Svm_sectioned', options, 3, 'Sectioned');

        cancerProb =  [0.69429197 0.66705237 0.39345565 1.         0.31506872 1.  0.34356235 0.66886696];
        predictedLabels = [1 1 0 1 0 1 0 1];
        getImageForClassifier(cancerProb, predictedLabels, coordinates, orderedRois, segmentMasks, sRGB, 'Knn_sectioned', options, 1, 'Sectioned');
        
        cancerProb = [0.52 0.56 0.52 0.58 0.61 0.53 0.51 0.44 ];
        predictedLabels = [1 1 1 1 1 1 1 0];
        getImageForClassifier(cancerProb, predictedLabels, coordinates, orderedRois, segmentMasks, sRGB, 'RF_sectioned', options, 2, 'Sectioned');
        
        trueLabels = [1 1 1 1 1 1 0 0];
        getImageForClassifier([0 0 0 0 0 0 0 0], trueLabels, coordinates, orderedRois, segmentMasks, sRGB, 'GD_sectioned', options, 2, 'Sectioned');

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
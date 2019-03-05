% Classify malignancy
validations = {'LeaveOut', 'Kfold'};
votingRules = {'majority', 'weighted majority', 'complex vote'};
frechet = @(Z1,ZJ) arrayfun(@(x) DiscreteFrechetDist(ZJ(x,:), Z1), (1:size(ZJ,1))');
distances = { frechet, 'correlation', 'euclidean'};
groups = {'fixed+unfixed', 'fixed', 'unfixed'}; %, 'unfixedleft' , 'goodleft', 'goodright'};
features = {'ref', 'refpca', 'reflda', 'refpcalda', 'lbp1', 'lbp2', 'reflbp1', 'reflbp2'}; %, 'spectrumlbp3'};
options.saveOptions.saveInHQ = true;
kernels = {'linear', 'rbf'};

if contains(lower(options.action), 'svm')
    classifier = 'svm'; 
    classifierSettings = kernels;
elseif contains(lower(options.action), 'knn')
    classifier = 'knn';
    classifierSettings = distances;
else 
     error('Unsupported classification method.Aborting...');
end
options.action = 'Classification';

dataset = 'estimated';

fprintf('Classifying %s data with %s classifier...\n', dataset, classifier);

m = 0;
validation = validations{2};
labelsAsText = true;

classifiers = struct('Input', {}, 'Features', {}, 'Validation', {}, 'VoteRule', {}, 'Neighbours', [], 'Setting', {}, 'FoldPerformance', [], 'Performance', [],...
    'Accuracy', [], 'AccuracySD', [], 'AUC', [], 'Fmeasure', []);
for g = 1:length(groups)
    [inputIdx, labels] = createClassifierInputIndexes(name, groups{g}, labelsAsText);
    [trainIdx, testIdx] = createCVTestIndexes(validation, labels);
    for p = 1:length(features)
        observations = classifierInput(dataset, inputIdx, labels, features{p}, name, groups{g});
        for i = 1:length(votingRules)
            for k = [1, 3, 5]
                for d = 1:length(classifierSettings)
                    [performance, cvPerformance] = crossValidation(trainIdx, testIdx, observations, labels, classifier, k, classifierSettings{d}, votingRules{i});
                    m = m + 1;
                    classifiers(m) = struct('Input', groups{g}, 'Features', features{p}, ...
                        'Validation', validation, 'VoteRule', votingRules{i}, 'Neighbours', k, 'Setting', classifierSettings{d}, 'FoldPerformance', cvPerformance,...
                        'Performance', performance, 'Accuracy', performance.Accuracy, 'AccuracySD', performance.AccuracySD, 'AUC', performance.AUC, 'Fmeasure', performance.Fmeasure);
                end
            end
        end
    end
end  
classifiers = orderfields(classifiers);
save( generateName(options, strcat(classifier,'_', dataset, '_', 'Classifier.mat')) , 'classifiers');

[~, sortIdx] = sort([classifiers.Accuracy], 'descend'); % or classificationError.AvgAUC
classifiers = classifiers(sortIdx);

if (options.showImages)
    %% Comparison plots 
    compareClassifiers(classifiers, options, classifier, dataset, validation)

    %% Comparison between input datasets 
    compareInputs(groups, classifiers, classifier, options, validation, dataset, 4);

    %% Comparison between features 
    compareFeatures(features, classifiers, classifier, options, validation, dataset, 5);
end

disp('Finished classification')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [performance] = getPerformanceStatistics(predictions, labels, scores, testIdxs)
%%Returns accuracy statistics for binary classification predictions
%Predictions: an 1D vector of prediction labels 
%Labels: an 1D vector of test labels 

    categories = {'Benign', 'Malignant'};
    if length(predictions) ~= length(labels) 
        error('Incompatible length of input data.')
    end
    
    if iscell(predictions) || iscell(labels)
        textLabels = labels; 
        textPredictions = predictions;
        labels = strcmp(labels, 'Malignant');
        predictions = strcmp(predictions, 'Malignant');
    else
        textLabels = categories(1 + labels);
        textPredictions = categories(1 + predictions);
    end
    % true = cancer condition is true OR the positive class is Malignant
	truePredictionIdxs = (predictions == labels);
	falsePredictionIdxs = (predictions ~= labels);
	
	a = sum(labels(truePredictionIdxs) == true);
	d = sum(labels(truePredictionIdxs) == false);
	c = sum(labels(falsePredictionIdxs) == true); 
	b = sum(labels(falsePredictionIdxs) == false);
	
	total = length(predictions);
	if (sum([a,b,c,d]) ~= total)
        error('Incorrect performance statistics.')
    end
	
    performance.TestIndexes = testIdxs;
    performance.Labels = textLabels;
    performance.Predictions = textPredictions;
	performance.ConfusionMatrix = [a, c ; b, d]; 
	performance.FalsePositiveRate = b / total * 100;
	performance.FalseNegativeRate = c / total * 100;
	performance.Sensitivity = a / (a + c) * 100; %or recall
	performance.Specificity = d / (b + d) * 100;
	performance.Precision = a / (a + b) * 100; %or positive predictive value
	performance.NegativePredictiveValue = d / (c + d) * 100;
	performance.Accuracy = (a + d) / total * 100;
    performance.Fmeasure = 2 * performance.Precision * performance.Sensitivity /(performance.Precision + performance.Sensitivity) / 100;
    [~,n] = size(scores);
    performance.Scores = scores(:,n);
    if any(labels == 0) && any(labels == 1)
        [performance.ROCX, performance.ROCY, performance.ROCT, performance.AUC] = perfcurve(textLabels, scores(:,n), 'Malignant'); % mark positive class
    else 
        performance.ROCX = [];
        performance.ROCY = [];
        performance.ROCT = [];
        performance.AUC = [];
    end
    performance = orderfields(performance);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [predictions, scores] = knn(testdata, database, databaselabels, k, distance, rule)
%%Knn classifies 'testdata' based on a 'database' and the respective 'databeselabels'
%testdata, database: rows = samples, columns = features
%Additional parameters: 
%Neighbors 'k'
%Distance 'distance'
%Voting rule 'rule'

    if (nargin < 4)
        k = 1;
    end
    if (nargin < 5)
        distance = 'euclidean';
    end
    if (nargin < 6)
        rule = 'majority';
    end
    
    [m, ~] = size(testdata);
    predictions = cell(m,1);
    scores = zeros(m,1);
    categories = {'Benign', 'Malignant'};
    for i = 1:m
        [knnIdx, D] = knnsearch(database, testdata(i,:), 'K', k, 'Distance', distance);
        malignantIdx = strcmp(databaselabels(knnIdx), 'Malignant'); % binary labels 
        nnrLabels = malignantIdx; % numerical labels as -1(benign) an +1 (malignant)
        nnrLabels(malignantIdx == 0) = -1;
        trueNeighbors = D(malignantIdx);
        falseNeighbors = D(~malignantIdx);
        
        if isempty(trueNeighbors) || isempty(falseNeighbors)
            rule = 'majority';
        end
        if strcmp(rule, 'complex vote')
            binaryLabel = -(4*log(min(falseNeighbors)) + log(mean(falseNeighbors)) + log(prod(falseNeighbors))) < ...
                -(4*log(min(trueNeighbors)) + log(mean(trueNeighbors)) + log(prod(trueNeighbors)));
            scores(i) = (4*log(min(falseNeighbors)) + log(mean(falseNeighbors)) + log(prod(falseNeighbors))) ...
                + (4*log(min(trueNeighbors)) + log(mean(trueNeighbors)) + log(prod(trueNeighbors))) * (-1);
            
        elseif strcmp(rule, 'weighted majority')
            binaryLabel = sum(1./falseNeighbors) < sum(1 ./ trueNeighbors);
            scores(i) =  sum(1./falseNeighbors) * (-1) + sum(1 ./ trueNeighbors) * 1;
            
        else %majority vote
            binaryLabel = round(mean(malignantIdx));
            scores(i) = mean(nnrLabels);
        end
        predictions(i) = categories(binaryLabel + 1);
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [predictions, performance] = classify(classifier, train, trainLabels, test, testLabels, k, setting, rule, testIndexes)
%%Classify returns the class predictions using 'classifier'
%test, train: rows = samples, columns = features

    if (nargin < 9)
        testIndexes = [];
    end
    
    if strcmp(classifier, 'knn')
        [predictions, scores] = knn(test, train, trainLabels, k, setting, rule); 
        
    elseif strcmp(classifier, 'svm')
        SVMModel = fitcsvm(train, trainLabels, 'KernelFunction', setting, 'Standardize',true,'ClassNames',{'Benign','Malignant'}, 'BoxConstraint', 1);
        [predictions,scores] = predict(SVMModel,test);
        %SVMModel = fitPosterior(SVMModel); [predictions,scores] = resubPredict(SVMModel);
        
    else 
       error('Unsupported classification method');
    end
    
    performance = getPerformanceStatistics(predictions, testLabels, scores, testIndexes);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [trainIdx, testIdx] = createCVTestIndexes(validation, labels)

    if strcmp(validation, 'Kfold')
        folds = 5;
        CVO = cvpartition(labels,'k',folds);
        
    elseif strcmp(validation, 'LeaveOut')
        CVO = cvpartition(labels, 'LeaveOut'); 
        
    else
        error('Unsupported validation method.');
    end
    
    folds = CVO.NumTestSets;
    testIdx = cell(folds, 1);
    trainIdx = cell(folds, 1);
    for i = 1:folds
        testIdx{i} = CVO.test(i);
        trainIdx{i} = CVO.training(i);
    end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [avgPerformance, cvPerformance] = crossValidation(trainIndexes, testIndexes, observations, labels, classifier, k, setting, rule)

%%CrossValidation returns crossvalidated accuracy statistics for
%%classification using 'classifier'
    [N, ~] = size(observations);  
    correctPredictions = zeros(N,1);
    totalScores = zeros(N,1);
    
    folds = length(trainIndexes);
    cvPerformance = struct('TestIndexes', [], 'Predictions', [], 'Labels', [],'Scores', [], 'ConfusionMatrix', [], ...
        'FalsePositiveRate', [], 'FalseNegativeRate', [], 'Sensitivity', [], 'Specificity', [], 'Precision', [], ...
        'NegativePredictiveValue', [], 'Accuracy', [], 'ROCX', [], 'ROCY', [], 'ROCT', [], 'AUC', [], 'Fmeasure', []);

    for i = 1:folds
        testIdx = testIndexes{i};
        trainIdx = trainIndexes{i};
        
        testLabels = labels(testIdx)';
        trainLabels = labels(trainIdx)';
        train = observations(trainIdx, :);
        test = observations(testIdx, :);  
        
        [predictions, performance] = classify(classifier, train, trainLabels, test, testLabels, k, setting, rule, find(testIdx == 1));
        cvPerformance(i) = performance;        
        correctPredictions(testIdx) = strcmp(predictions, testLabels);
        totalScores(testIdx) = cvPerformance(i).Scores;
    end
    
    %% avarage data 
    avgPerformance.Labels = {labels};
    avgPerformance.Scores = totalScores;
    avgPerformance.IsTruePrediction = correctPredictions;
    avgPerformance.FalsePositiveRate =  mean([cvPerformance.FalsePositiveRate]);
    avgPerformance.FalseNegativeRate =  mean([cvPerformance.FalseNegativeRate]);
    avgPerformance.Sensitivity =  mean([cvPerformance.Sensitivity]);
    avgPerformance.Specificity =  mean([cvPerformance.Specificity]);
    avgPerformance.Precision =  mean([cvPerformance.Precision]);
    avgPerformance.Accuracy =  mean([cvPerformance.Accuracy]);
    avgPerformance.AccuracySD =  std([cvPerformance.Accuracy]);
    avgPerformance.Fmeasure = mean([cvPerformance.Fmeasure]);
    
    nonEmptyIdx = ~cellfun('isempty', {cvPerformance.ROCY});
    if sum(cellfun('isempty', {cvPerformance.ROCY})) > 0
        emptyClassifiers = sum(cellfun('isempty', {cvPerformance.ROCY}));
        fprintf('Empty classifiers: %d\n', emptyClassifiers);
    end
    [avgPerformance.ROCX, avgPerformance.ROCY, avgPerformance.ROCT, AUC] = perfcurve({cvPerformance(nonEmptyIdx).Labels}, {cvPerformance(nonEmptyIdx).Scores}, 'Malignant', 'XVals', 'all');
    avgPerformance.AUC = AUC(1);

        
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [inputIdx, labels] = createClassifierInputIndexes(name, criterion, labelsAsText)

    if (nargin < 3)
        labelsAsText = false;
    end
    
    load( fullfile(generateName([], 'input'), name, 'ID.mat'), 'ID');
    
    if strcmp(criterion, 'unique') || strcmp(criterion, 'fixed+unfixed')
        [~, inputIdx, ~] = unique(strcat({ID.SpectrumFile}, {ID.T}), 'last');
        
    elseif strcmp(criterion, 'fixed+unfixed')  % not cut 
        [~, unIdx, ~] = unique(strcat({ID.SpectrumFile}, {ID.T}), 'last');
        inputIdx = intersect(unIdx, find([ID.IsCut] == false));
        
    elseif strcmp(criterion, 'fixed')
        [~, unIdx, ~] = unique(strcat({ID.SpectrumFile}, {ID.T}), 'last');
        inputIdx = intersect(unIdx, find([ID.IsFixed] == true));

    elseif strcmp(criterion, 'unfixed') ||  strcmp(criterion, 'unfixedright') 
        [~, unIdx, ~] = unique(strcat({ID.SpectrumFile}, {ID.T}), 'last');
        inputIdx = intersect(unIdx, find([ID.IsFixed] == false));
        
%     elseif strcmp(criterion, 'unfixedleft')
%         [~, unIdx, ~] = unique(strcat({ID.SpectrumFile}, {ID.T}), 'first');
%         inputIdx = intersect(unIdx, find([ID.IsFixed] == false));
        
%     elseif strcmp(criterion, 'goodright') || strcmp(criterion, 'good')
%         [~, unIdx, ~] = unique(strcat({ID.SpectrumFile}, {ID.T}), 'last');
%         goodIdx = intersect(unIdx , find([ID.IsGood]));
%         %goodIdx = union(find(strcmp([ID.Sample], '9913')), union(find(strcmp([ID.Sample], '9933')), union(find(strcmp([ID.Sample], '9940')), find(strcmp([ID.Sample], '9956')))));
%         inputIdx = intersect(unIdx, goodIdx);
%     
%     elseif strcmp(criterion, 'goodleft')
%         [~, unIdx, ~] = unique(strcat({ID.SpectrumFile}, {ID.T}), 'first');
%         goodIdx = intersect(unIdx , find([ID.IsGood]));
%         %goodIdx = union(find(strcmp([ID.Sample], '9913')), union(find(strcmp([ID.Sample], '9933')), union(find(strcmp([ID.Sample], '9940')), find(strcmp([ID.Sample], '9956')))));
%         inputIdx = intersect(unIdx, goodIdx);
        
    elseif strcmp(criterion, 'all')
        inputIdx = 1:numel(ID);      
        
    else
        error('Not implemented yet.')
    end
    
    if ~ strcmp(criterion, 'all')
        load(fullfile(generateName([], 'output'), name, 'ReflectanceEstimationPreset', 'out.mat'), 'goodSampleIndexes');
        inputIdx = intersect(inputIdx, find(goodSampleIndexes));
    end
        
    A = [ID.IsBenign]';
    if (labelsAsText)
        X = {'Malignant', 'Benign'};
        labels = X(1+A);
    else
        labels = ~A;
    end
    labels = labels(inputIdx);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function observations = classifierInput(dataset, inputIdx, labels, features, name, inputgroup)

    if (nargin < 6)
        inputgroup = '';
    end

    reflectanceFeats = [];
    if contains(features, 'ref')
        if strcmp(dataset, 'measured')
            load(fullfile(generateName([], 'input'), name, 'in.mat'), 'Spectra');
            Gall = Spectra;

        elseif strcmp(dataset, 'estimated')
            load(fullfile(generateName([], 'output'), name, 'ReflectanceEstimationPreset', 'out.mat'), 'EstimatedSpectra');
            Gall = EstimatedSpectra;

        else
            error('Not acceptable input. Choose "measured" or "estimated".')
        end

        reflectanceFeats = Gall(inputIdx, :);  % G rows are observations and columns are variables

        if contains(features, 'pca' ) || contains(features, 'lda')
            parts = split(features, 'lbp');
            projection = strrep(parts(1), 'ref', '');
            [~, scores] = dimensionReduction(projection, reflectanceFeats, labels);
            reflectanceFeats = scores(:, 1:10);
        end
    end
    
    lbpFeats = [];
    if contains(features, 'lbp')
        load( fullfile(generateName([], 'output'), name, 'ReflectanceEstimationPreset', 'out.mat'), 'MultiScaleLbpFeatures');
        if contains(features, 'lbp1')
            scales = 1;
        elseif contains(features, 'lbp2')
            scales = 2;
        elseif contains(features, 'lbp3')
            scales = 3;
        else
            scales = length(MultiScaleLbpFeatures);
        end      
        lbps = size(MultiScaleLbpFeatures{1},2);
        lbpFeats = zeros(length(labels), scales * lbps);      
        for i = 1:scales
            lbpFeatures = MultiScaleLbpFeatures{i};
            rangeIdx = (i-1)*lbps + (1:lbps);
            lbpFeats(:,rangeIdx) = lbpFeatures(inputIdx,:);
        end
    end    
    observations = [reflectanceFeats, lbpFeats]; 
    
    %% Export feature set 
    featuresFilename = fullfile(generateName([], 'output'), name, 'Classification', strcat(strjoin({'FeatureVectors', inputgroup, features }, '_'), '.mat' ));
    if ~exist( featuresFilename, 'file')
        save(featuresFilename, 'observations', 'labels');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = plotROC(selectedClassifier, options, classifier, dataset, fig, name)
    figTitle = strcat('ROC for [', selectedClassifier.Input,'+', selectedClassifier.Features ,'] dataset.',...
        'Accuracy:', num2str(selectedClassifier.Accuracy), '\pm', num2str(selectedClassifier.AccuracySD));
    options.saveOptions.plotName = generateName(options, strjoin({classifier, dataset, name}, '_')); 
    plots('roc', fig, [], '', 'Performance', selectedClassifier.Performance , 'FoldPerformance', selectedClassifier.FoldPerformance, ...
        'SaveOptions', options.saveOptions, 'Title', figTitle);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = compareClassifiers(classifiers, options, classifier, dataset, validation)

    [~, maxAccurIdx] = max([classifiers.Accuracy]);
    classifiers(maxAccurIdx)

    if ~contains(validation, 'LeaveOut')
        [~, maxAUCIdx] = max([classifiers.AUC]);
        classifiers(maxAUCIdx)

        [~, maxF1Idx] = max([classifiers.Fmeasure]);
        classifiers(maxF1Idx)

        plotROC( classifiers(maxAccurIdx), options, classifier, dataset, 1, 'MaxAccuracyRoc');
        plotROC( classifiers(maxAUCIdx), options, classifier, dataset, 2, 'MaxAUCRoc');
        plotROC( classifiers(maxF1Idx), options, classifier, dataset, 3, 'MaxF1Roc');
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = compareInputs(groups, classifiers, classifier, options, validation, dataset, fig)

    accur = zeros(1,length(groups));
    auc = zeros(1,length(groups));
    for i = 1:length(groups)
        classIdx = strcmp({classifiers.Input}, groups{i});
        selectedClassifiers = classifiers(classIdx);
        if contains(validation, 'LeaveOut')
            [~, maxIdx] = max([selectedClassifiers.Accuracy]);
            selectedClassifier = selectedClassifiers(maxIdx)
            auc(i) = NaN;
            accur(i) = selectedClassifier.Accuracy;
        else
            [~, maxIdx] = max([selectedClassifiers.AUC]);
            selectedClassifier = selectedClassifiers(maxIdx)
            auc(i) = selectedClassifier.AUC;
            accur(i) = selectedClassifier.Accuracy;
        end
   
    end
    
    figTitle = sprintf('%s Classifier Performance for Different Input Datasets', upper(classifier));
    options.saveOptions.plotName = generateName(options, strjoin({classifier, dataset, validation, 'compareInput'}, '_'));
    plots('performanceComparison', fig, [], '', 'LineNames', groups, 'Performance', [auc; accur], 'SaveOptions', options.saveOptions, 'Title', figTitle);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = compareFeatures(features, classifiers, classifier, options, validation, dataset, fig)

    accur = zeros(1,length(features));
    auc = zeros(1,length(features));
    for i = 1:length(features)
        classIdx = strcmp({classifiers.Features}, features{i});
        selectedClassifiers = classifiers(classIdx);
        if contains(validation, 'LeaveOut')
            [~, maxIdx] = max([selectedClassifiers.Accuracy]);
            selectedClassifier = selectedClassifiers(maxIdx)
            auc(i) = NaN;
            accur(i) = selectedClassifier.Accuracy;
        else
            [~, maxIdx] = max([selectedClassifiers.AUC]);
            selectedClassifier = selectedClassifiers(maxIdx)
            auc(i) = selectedClassifier.AUC;
            accur(i) = selectedClassifier.Accuracy;
        end
   
    end
    
    figTitle = sprintf('%s Classifier Performance for Different Feature Sets', upper(classifier));
    options.saveOptions.plotName = generateName(options, strjoin({classifier, dataset, validation, 'compareFeatures'}, '_'));
    plots('performanceComparison', fig, [], '', 'LineNames', features, 'Performance', [auc; accur], 'SaveOptions', options.saveOptions, 'Title', figTitle);
    
end
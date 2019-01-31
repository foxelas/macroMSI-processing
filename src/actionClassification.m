% Classify malignancy
if contains(lower(options.action), 'svm')
    classifier = 'svm'; 
elseif contains(lower(options.action), 'knn')
    classifier = 'knn';
else 
     error('Unsupported classification method.Aborting...');
end

version = 'measured';

fprintf('Classifying %s data with %s classifier...\n', version, classifier);

validations = {'LeaveMOut', 'Kfold'};
votingRules = {'majority', 'weighted majority', 'complex vote'};
frechet = @(Z1,ZJ) arrayfun(@(x) DiscreteFrechetDist(ZJ(x,:), Z1), (1:size(ZJ,1))');
distances = { frechet, 'correlation', 'chebychev', 'euclidean'};
groups = {'unique', 'fixed', 'unfixed'};
projections = {'spectrum'};%, 'pca', 'lda', 'pcalda', 'spectrumlbp'};
options.saveOptions.saveInHQ = true;
kernels = {'linear', 'rbf'};

m = 0;
validation = validations{2};
labelsAsText = true;


%% KNN
if strcmp(classifier, 'knn')   
    classificationError = struct('Input', {}, 'Projection', {}, 'Validation', {}, 'VoteRule', {}, 'Neighbours', [], 'Distance', {}, 'FoldPerformance', [], 'Performance', [],...
        'Accuracy', [], 'AccuracySD', [], 'AUC', []);
    for g = 1:length(groups)
        for p = 1:length(projections)
            [Gun, labels] = classifierInput(version, groups{g}, projections{p}, name, labelsAsText);
            for i = 1:length(votingRules)
                for k = [1, 3, 5]
                    for d = 1:length(distances)

                        [performance, cvPerformance] = crossValidation(validation, Gun, labels, classifier, k,  distances{d}, votingRules{i}, []);                          
                        m = m + 1;
                        classificationError(m) = struct('Input', groups{g}, 'Projection', projections{p}, ...
                            'Validation', validation, 'VoteRule', votingRules{i}, 'Neighbours', k, 'Distance', distances{d}, 'FoldPerformance', cvPerformance,...
                            'Performance', performance, 'Accuracy', performance.Accuracy, 'AccuracySD', performance.AccuracySD, 'AUC', performance.AUC);
                    end
                end
            end
        end
    end  
    
%% SVM
elseif strcmp(classifier, 'svm')

    classificationError = struct('GenLoss', [], 'Input', {}, 'Projection', {}, 'Validation', {}, 'Kernel', {}, 'FoldPerformance', [], 'Performance', [],...
        'Accuracy', [], 'AccuracySD', [], 'AUC', []);
    for g = 1:length(groups)
        for p = 1:length(projections)
            [Gun, labels] = classifierInput(version, groups{g}, projections{p}, name, labelsAsText);
            for i = 1:length(kernels)

                [performance,cvPerformance] = crossValidation(validation, Gun, labels, classifier, [], [], [], kernels{i});                                             
                m = m + 1;
                classificationError(m) = struct('Input', groups{g}, 'Projection', projections{p},  ...
                    'Validation', validation, 'Kernel', kernels{i}, 'FoldPerformance', cvPerformance, 'GenLoss', [], ...
                    'Performance', performance, 'Accuracy', performance.Accuracy, 'AccuracySD', performance.AccuracySD, 'AUC', performance.AUC);
            end
        end
    end
        
end

classificationError = orderfields(classificationError);
dirr = fullfile(options.saveOptions.savedir, 'Classification');
if ~exist(dirr, 'dir')
    mkdir(dirr);
    addpath(dirr);
end
save( fullfile(dirr, strcat(classifier,'_', version, '_', 'Error.mat')),'classificationError');

[~, sortIdx] = sort([classificationError.Accuracy], 'descend'); % or classificationError.AvgAUC
classificationError = classificationError(sortIdx);

[~, maxAccurIdx] = max([classificationError.Accuracy]);
figTitle = strcat('ROC for [', classificationError(maxAccurIdx).Input,'+', classificationError(maxAccurIdx).Projection ,'] dataset.',...
    'Accuracy:', num2str(classificationError(maxAccurIdx).Accuracy), '\pm', num2str(classificationError(maxAccurIdx).AccuracySD));
options.saveOptions.plotName = fullfile( options.saveOptions.savedir, 'Classification', strcat(classifier,'_', version, '_MaxAccuracyRoc'));
plots('roc', 1, [], '', 'Performance', classificationError(maxAccurIdx).Performance , 'FoldPerformance', classificationError(maxAccurIdx).FoldPerformance, 'SaveOptions', options.saveOptions, 'Title', figTitle);

[~, maxAUCIdx] = max([classificationError.AUC]);
figTitle = strcat('ROC for [', classificationError(maxAUCIdx).Input,'+', classificationError(maxAUCIdx).Projection ,'] dataset.',...
    'Accuracy:', num2str(classificationError(maxAUCIdx).Accuracy), '\pm', num2str(classificationError(maxAUCIdx).AccuracySD));
options.saveOptions.plotName = fullfile( options.saveOptions.savedir, 'Classification', strcat(classifier, '_', version, '_MaxAUCRoc'));
plots('roc', 2, [], '', 'Performance', classificationError(maxAUCIdx).Performance , 'FoldPerformance', classificationError(maxAUCIdx).FoldPerformance, 'SaveOptions', options.saveOptions, 'Title', figTitle);
%options.saveOptions.plotName = generateName(options, ['Classification error of ', version, ' spectra with ', validation]);
%plots('classificationErrors', 2, [], '', 'Errors', classificationError, 'SaveOptions', options.saveOptions)
    

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
		error('Incorrect performance statistics.');
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

function [predictions, scores, performance] = classify(classifier, train, trainLabels, test, testLabels, k, distance, rule, kernel, testIndexes)
%%Classify returns the class predictions using 'classifier'
%test, train: rows = samples, columns = features

    if (nargin < 10)
        testIndexes = [];
    end
    
    if strcmp(classifier, 'knn')
        [predictions, scores] = knn(test, train, trainLabels, k, distance, rule); 
        
    elseif strcmp(classifier, 'svm')
        SVMModel = fitcsvm(train, trainLabels, 'KernelFunction', kernel, 'Standardize',true,'ClassNames',{'Benign','Malignant'}, 'BoxConstraint', 1);
        [predictions,scores] = predict(SVMModel,test);
        %SVMModel = fitPosterior(SVMModel); [predictions,scores] = resubPredict(SVMModel);
        
    else 
       error('Unsupported classification method');
    end
    
    performance = getPerformanceStatistics(predictions, testLabels, scores, testIndexes);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [avgPerformance, cvPerformance] = crossValidation(validation, data, labels, classifier, k, distance, rule, kernel)
%%CrossValidation returns crossvalidated accuracy statistics for
%%classification using 'classifier'
    if (nargin < 3)
        validation = 'Kfold';
    end
        
    [m, ~] = size(data);
    if strcmp(validation, 'Kfold')
        folds = 10;
        indices = crossvalind('Kfold', m, folds);
        
    elseif strcmp(validation, 'LeaveMOut')
        folds = m;
        indices = crossvalind('LeaveMOut', m, 1); 
        
    else
        error('Unsupported validation method.');
    end
    
    correctPredictions = zeros(m,1);
    totalScores = zeros(m,1);
    
    cvPerformance = struct('TestIndexes', [], 'Predictions', [], 'Labels', [],'Scores', [], 'ConfusionMatrix', [], ...
        'FalsePositiveRate', [], 'FalseNegativeRate', [], 'Sensitivity', [], 'Specificity', [], ...
        'Precision', [], 'NegativePredictiveValue', [], 'Accuracy', [], 'ROCX', [], 'ROCY', [], 'ROCT', [], 'AUC', []);
    
    for i = 1:folds
        
        if strcmp(validation, 'Kfold')
            testIdx = (indices == i);
            trainIdx = ~testIdx;
            
        elseif strcmp(validation, 'LeaveMOut')
            trainIdx = true(size(labels));
            testIdx = false(size(labels));
            trainIdx(m) = false;
            testIdx(m) = true;
        end
        
        testLabels = labels(testIdx);
        trainLabels = labels(trainIdx);
        train = data(trainIdx, :);
        test = data(testIdx, :);
        
        [predictions, scores, performance] = classify(classifier, train, trainLabels, test, testLabels, k, distance, rule, kernel, find(testIdx == 1));
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
    nonEmptyIdx = ~cellfun('isempty', {cvPerformance.ROCY});
    [avgPerformance.ROCX, avgPerformance.ROCY, avgPerformance.ROCT, auc] = perfcurve({cvPerformance(nonEmptyIdx).Labels}, {cvPerformance(nonEmptyIdx).Scores}, 'Malignant', 'XVals', 'all');
    avgPerformance.AUC = auc(1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [G, labels] = classifierInput(version, group, features, name, labelsAsText)

    if (nargin < 5)
        labelsAsText = false;
    end

    [Gfx, ~, subIdx, labelsfx] = subset(version, name, group);

    if contains(features, 'pca' ) || contains(features, 'lda')
        [~, scores] = dimensionReduction(features, Gfx, double(labelsfx));
        G = scores(:, 1:10);

    else %contains(features, 'spectrum')
        G = Gfx;
    end

    if contains(features, 'lbp')
        e =  matfile(fullfile('..', '..', 'output', name, 'reflectanceestimationpreset', 'out.mat'));
        multiScaleLbpFeatures = e.multiScaleLbpFeatures;
        scales = length(multiScaleLbpFeatures);
        lbps = size(multiScaleLbpFeatures{1},2);
        Gplus = zeros(length(labelsfx), length(G) + scales * lbps);
        Gplus(:,1:size(G, 2)) = G;
        for i = 1:scales
            lbpFeatures = multiScaleLbpFeatures{i};
            rangeIdx = length(G) + (i-1)*lbps + (1:lbps);
            Gplus(:,rangeIdx) = lbpFeatures(subIdx,:);
        end
        G = Gplus;
    end

    labels = ~labelsfx';

    if (labelsAsText)
        X = {'Benign', 'Malignant'};
        labels = X(1+labels)';
    end

end

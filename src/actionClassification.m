% Classify malignancy
if contains(lower(options.action), 'svm')
    classifier = 'svm'; 
elseif contains(lower(options.action), 'knn')
    classifier = 'knn';
else 
     error('Unsupported classification method.Aborting...');
end

version = 'estimated';

fprintf('Classifying %s data with %s classifier...\n', version, classifier);

validations = {'LeaveMOut', 'Kfold'};
votingRules = {'majority', 'weighted majority', 'complex vote'};
distances = {'correlation', 'chebychev', 'euclidean'};
groups = {'unique', 'fixed', 'unfixed', 'good'};
projections = {'spectrum'}; %, 'pca', 'lda', 'pcalda', 'spectrumlbp'};
options.saveOptions.saveInHQ = true;
kernels = {'linear', 'rbf'};

m = 0;
validation = validations{2};

%% KNN
if strcmp(classifier, 'knn')   
    classificationError = struct('Accuracy', [], 'FalsePositives', [], 'FalseNegatives', [], 'Input', {}, 'Projection', {}, 'Validation', {}, 'VoteRule', {}, 'Neighbours', [], 'Distance', {});
    for g = 1:length(groups)
        for p = 1:length(projections)
            [Gun, labels] = classifierInput(version, groups{g}, projections{p}, name);
            for i = 1:length(votingRules)
                for k = [1, 3, 5]
                    for d = 1:length(distances)

                        [a, b, c] = crossValidation(validation, Gun, labels, classifier, k,  distances{d}, votingRules{i}, []);                          
                        m = m + 1;
                        classificationError(m) = struct('Accuracy', a, 'FalsePositives', b, 'FalseNegatives', c, 'Input', groups{g}, 'Projection', projections{p}, ...
                            'Validation', validation, 'VoteRule', votingRules{i}, 'Neighbours', k, 'Distance', distances{d});
                    end
                end
            end
        end
    end  
    
%% SVM
elseif strcmp(classifier, 'svm')

    labelsAsText = true;
    classificationError = struct('Accuracy', [], 'FalsePositives', [], 'FalseNegatives', [], 'GenLoss', [], 'Input', {}, 'Projection', {}, 'Validation', {}, 'Kernel', {});
    for g = 1:length(groups)
        for p = 1:length(projections)
            [Gun, labels] = classifierInput(version, groups{g}, projections{p}, name, labelsAsText);
            for i = 1:length(kernels)

                [a, b, c] = crossValidation(validation, Gun, labels, classifier, [], [], [], kernels{i});                                             
                m = m + 1;
                classificationError(m) = struct('Accuracy', a, 'FalsePositives', b, 'FalseNegatives', c, 'Input', groups{g}, 'Projection', projections{p}, ...
                            'Validation', validation, 'Kernel', kernels{i});
            end
        end
    end
        
end
 
save('classificationError.mat', 'classificationError');

pause;

classificationErrorFields = fieldnames(classificationError);
classificationErrorCell = struct2cell(classificationError);
sz = size(classificationErrorCell);
% Convert to a matrix
classificationErrorCell = reshape(classificationErrorCell, sz(1), []);
% Make each field a column
classificationErrorCell = classificationErrorCell';
% Sort by first field "name"
classificationErrorCell = sortrows(classificationErrorCell, -1);
classificationErrorCell = reshape(classificationErrorCell', sz);
classificationError = cell2struct(classificationErrorCell, classificationErrorFields, 1);
options.saveOptions.plotName = generateName(options, ['Classification error of ', version, ' spectra with ', validation]);
plots('classificationErrors', 2, [], '', 'errors', classificationError, 'saveOptions', options.saveOptions)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [accuracyRate, falsePositiveRate, falseNegativeRate] = getAccuracyStatistics(predictions, labels)
%%Returns accuracy statistics for binary classification predictions
%Predictions: an 1D vector of prediction labels 
%Labels: an 1D vector of test labels 

    if length(predictions) ~= length(labels) 
        error('Incompatible length of input data.')
    end
    % true = cancer condition is true
    falsePredictionIdxs = ~(predictions == labels);
    falsePositives = sum(labels(falsePredictionIdxs) == false);
    falseNegatives = sum(labels(falsePredictionIdxs) == true);
    
    total = length(predictions);
    falsePositiveRate = falsePositives / total * 100;
    falseNegativeRate = falseNegatives / total * 100;
    accuracyRate = 100 - falsePositiveRate - falseNegativeRate;
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
    
    [m, n] = size(testdata);
    predictions = zero(m,1);
    scores = zero(m,1);
    for i = 1:m
        [knnIdx, D] = knnsearch(testdata(i,:), database, 'K', k, 'Distance', distance);
        nnrLabels = databaselabels(knnIdx);
        trueNeighbors = D(nnrLabels);
        falseNeighbors = D(~nnrLabels);

        if sum(nnrLabels) == 0
            predictions(i) = 0;
            scores(i) = -k;
            
        elseif sum(nnrLabels) == k
            predictions(i) = 1;
            scores(i) = k;
            
        elseif strcmp(rule, 'complex vote')
            predictions(i) = (4*log(min(falseNeighbors)) + log(mean(falseNeighbors)) + log(prob(falseNeighbors))) > ...
                (4*log(min(trueNeighbors)) + log(mean(trueNeighbors)) + log(prob(trueNeighbors)));
            scores(i) = (4*log(min(falseNeighbors)) + log(mean(falseNeighbors)) + log(prob(falseNeighbors))) - ...
                (4*log(min(trueNeighbors)) + log(mean(trueNeighbors)) + log(prob(trueNeighbors)));
            
        elseif strcmp(rule, 'weighted majority')
            predictions(i) = sum(1./falseNeighbors) < sum(1 ./ trueNeighbors);
            scores(i) = sum(1./falseNeighbors) - sum(1 ./ trueNeighbors);
            
        else %majority vote
            predictions(i) = round(mean(nnrLabels));
            scores(i) = trueNeighbors - falseNeighbors; 
        end
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [predictions, scores, accuracyRate, falsePositiveRate, falseNegativeRate] = classify(classifier, train, trainLabels, test, testLabels, k, distance, rule, kernel)
%%Classify returns the class predictions using 'classifier'
%test, train: rows = samples, columns = features

    if strcmp(classifier, 'knn')
        [predictions, scores] = knn(test, train, trainLabels, k, distance, rule); 
        
    elseif strcmp(classifier, 'svm')
        SVMModel = fitcsvm(train, trainLabels, 'KernelFunction', kernel, 'Standardize',true,'ClassNames',{'Benign','Malignant'}, 'BoxConstraint', 1);
        [predictions,scores] = predict(SVMModel,test);
        
    else 
       error('Unsupported classification method');
    end
    
    [accuracyRate, falsePositiveRate, falseNegativeRate] = getAccuracyStatistics(predictions, testLabels);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [avgAccuracyRate, avgFalsePositiveRate, avgFalseNegativeRate] = crossValidation(validation, data, labels, classifier, k, distance, rule, kernel)
%%CrossValidation returns crossvalidated accuracy statistics for
%%classification using 'classifier'
    if (nargin < 3)
        validation = 'Kfold';
    end
    
    avgAccuracyRate = 0;
    avgFalsePositiveRate = 0;
    avgFalseNegativeRate = 0;
        
    [m, n] = size(testdata);
    if strcmp(validation, 'Kfold')
        folds = 10;
        indices = crossvalind('Kfold', n, folds);
        
    elseif strcmp(validation, 'LeaveMOut')
        folds = n;
        indices = crossvalind('LeaveMOut', n, 1); 
        
    else
        error('Unsupported validation method.');
    end
        
    for i = 1:folds
        
        if strcmp(validation, 'Kfold')
            testIdx = (indices == i);
            trainIdx = ~testIdx;
            
        elseif strcmp(validation, 'LeaveMOut')
            trainIdx = true(size(labels));
            testIdx = false(size(labels));
            trainIdx(n) = false;
            testIdx(n) = true;
        end
        testLabels = labels(testIdx);
        trainLabels = labels(trainIdx);
        train = data(trainIdx, :);
        test = data(testIdx, :);
        
        [predictions, scores, accuracyRate, falsePositiveRate, falseNegativeRate] = classify(classifier, train, trainLabels, test, testLabels, k, distance, rule, kernel);
        avgAccuracyRate = avgAccuracyRate + accuracyRate / folds;
        avgFalseNegativeRate = avgFalseNegativeRate + falseNegativeRate / folds;
        avgFalsePositiveRate = avgFalsePositiveRate + falsePositiveRate / folds;
    end
    
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
        e = matfile(fullfile('..', '..', 'output', name, 'LBP', 'lbp.mat'));
        lbpFeatures = e.lbpFeatures;
        G = [G lbpFeatures(subIdx,:)];
    end

    labels = ~labelsfx';

    if (labelsAsText)
        X = {'Benign', 'Malignant'};
        labels = X(1+labels);
    end

end

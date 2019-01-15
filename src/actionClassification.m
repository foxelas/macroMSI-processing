% Classify malignancy
if contains(action, 'svm')
    method = 'svm'; 
elseif contains(action, 'knn')
    method = 'knn';
elseif contains(action, 'lbp')
    method = 'lbp';
else 
     error('Unsupported classification method.Aborting...');
end

version = 'estimated';

fprintf('Classifying %s data with %s classifier...\n', version, method);

validations = {'LeaveMOut', 'Kfold'};
validationShort = {'1O', 'Kf'};
votingRules = {'majority', 'weighted majority', 'complex vote'};
votingRulesAbbr = {'M', 'W', 'C'};
distances = {'correlation', 'chebychev', 'euclidean'};
distancesAbbr = {'Cr', 'Ch', 'Eu'};
groups = {'unique', 'fixed', 'unfixed', 'good'};
projections = {'spectrum'}%, 'pca', 'lda', 'pcalda'};
options.saveOptions.saveInHQ = true;
kernels = {'linear', 'rbf'};


if strcmp(method, 'knn')

    for j = 1:length(validations)
        m = 0;
        %% KNN
        classificationError = struct('Accuracy', [], 'TypeI', [], 'TypeII', [], 'Input', {}, 'Projection', {}, 'Validation', {}, 'VoteRule', {}, 'Neighbours', [], 'Distance', {});
        for g = 1:length(groups)
            for p = 1:length(projections)
                [Gun, labels] = classifierInput(version, groups{g}, projections{p}, name);
                for i = 1:length(votingRules)
                    for k = [1, 3, 5]
                        for d = 1:length(distances)
                            
                            [a, b, c] = knn(Gun, labels, k, distances{d}, votingRules{i}, validations{j}); % labels = 1 for cancer
                            
                            m = m + 1;
                            classificationError(m) = struct('Accuracy', a, 'TypeI', b, 'TypeII', c, 'Input', groups{g}, 'Projection', projections{p}, ...
                                'Validation', validations{j}, 'VoteRule', votingRules{i}, 'Neighbours', k, 'Distance', distances{d});
                        end
                    end
                end
            end
        end
        
    end
    
        
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
    options.saveOptions.plotName = generateName(options, ['Classification error of ', version, ' spectra with ', validations{j}]);
    plots('classificationErrors', 2, [], '', 'errors', classificationError, 'saveOptions', options.saveOptions)

elseif strcmp(method, 'svm')

    labelsAsText = true;
    for j = 2:length(validations)
        m = 0;
        %% SVM
        classificationError = struct('Accuracy', [], 'TypeI', [], 'TypeII', [], 'GenLoss', [], 'Input', {}, 'Projection', {}, 'Validation', {}, 'Kernel', {});
        for g = 1:length(groups)
            for p = 1:length(projections)
                [Gun, labels] = classifierInput(version, groups{g}, projections{p}, name, labelsAsText);
                for i = 1:length(kernels)
                    
                    [a, b, c] = getSVMAccuracy(Gun, labels, kernels{i}, validations{j});
                    
                    SVMModel = fitcsvm(Gun, labels, 'KernelFunction', kernels{j}, 'Standardize',false,'ClassNames',{'Benign','Malignant'});
                    if (j ==1)
                        %leave out
                        CVSVMModel = crossval(SVMModel,'Holdout',0.15);
                        genloss = kfoldLoss(CVSVMModel);
                    elseif (j == 2) 
                        %10fold cross
                        CVSVMModel = crossval(SVMModel);
                        genloss = kfoldLoss(CVSVMModel);
                    end
                   
                    m = m + 1;
                    classificationError(m) = struct('Accuracy', a, 'TypeI', b, 'TypeII', c, 'GenLoss', genloss, 'Input', groups{g}, 'Projection', projections{p}, ...
                                'Validation', validations{j}, 'Kernel', kernels{i});
                end
            end
        end
    end
elseif strcmp(method, 'lbp')
    
    load('lbp.mat');
    classificationError = struct('Accuracy', [], 'TypeI', [], 'TypeII', [], 'GenLoss', [], 'Validation', {}, 'Kernel', {});

    j = 2;
    m = 0;
    for i = 1:length(kernels)
        [a, b, c] = getSVMAccuracy(lbpFeatures, lbpLabels, kernels{i}, validations{j});
        SVMModel = fitcsvm(lbpFeatures, lbpLabels, 'KernelFunction', kernels{j}, 'Standardize',false,'ClassNames',{'Benign','Malignant'});
        CVSVMModel = crossval(SVMModel);
        genloss = kfoldLoss(CVSVMModel);
        
        m = m + 1;
        classificationError(m) = struct('Accuracy', a, 'TypeI', b, 'TypeII', c, 'GenLoss', genloss, ...
            'Validation', validations{j}, 'Kernel', kernels{i});
    end
 
    save('classificationError.mat', 'classificationError');
end
    
function [a, falsePositivesRate, falseNegativesRate] = getSVMAccuracy(data, labels, kernel, validation)
                    
if strcmp(validation, 'LeaveMOut')
    r = zeros(1, 2);
    for n = 1:length(labels)
        trainIdx = true(size(labels));
        testIdx = false(size(labels));
        trainIdx(n) = false;
        testIdx(n) = true;
        testLabel = labels(testIdx);
        trainLabels = labels(trainIdx);
        train = data(trainIdx, :);
        test = data(testIdx, :);
        
        SVMModel = fitcsvm(train, trainLabels, 'KernelFunction', kernel, 'Standardize', true, 'ClassNames',{'Benign','Malignant'}, 'BoxConstraint', 1);
        [prediction,score] = predict(SVMModel,test);
        
        r(1) = r(1) + double(strcmp(prediction, testLabel));
        r(2) = r(2) + double(~strcmp(prediction, testLabel) && strcmp(testLabel, 'Malignant')); % type I error       
    end
    r = r ./ length(labels);
    
elseif strcmp(validation, 'Kfold')
    indices = crossvalind('Kfold', length(labels), 10);
    r = zeros(1, 3);
    for i = 1:10
        testIdx = (indices == i);
        trainIdx = ~testIdx;
        testLabels = labels(testIdx);
        trainLabels = labels(trainIdx);
        train = data(trainIdx, :);
        test = data(testIdx, :);
        rf = zeros(1, 3);
        for j = 1:sum(testIdx)
            SVMModel = fitcsvm(train, trainLabels, 'KernelFunction', kernel, 'Standardize',true,'ClassNames',{'Benign','Malignant'}, 'BoxConstraint', 1);
            [prediction,score] = predict(SVMModel,test(j,:));
            
            rf(1) = rf(1) + double(strcmp(prediction, testLabels(j)));
            rf(2) = rf(2) + double(~strcmp(prediction, testLabels(j)) && strcmp(testLabels(j), 'Malignant')); % type I error
        end
        r = r + rf ./ (sum(testIdx) * 10);
        
    end
    
else
    error('Not implemented')
end
r = r * 100;
[a, falsePositivesRate, falseNegativesRate] = deal(r(1), r(2), 100 - sum(r));

end
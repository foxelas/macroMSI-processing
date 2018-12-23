function [a, falsePositivesRate, falseNegativesRate] = knn( data, labels, k, distance, rule, validation)
% Label tags: 0 is normal , 1 is cancer 

if (nargin < 3)
    k = 1;
end
if (nargin < 4)
    distance = 'euclidean';
end
if (nargin < 5)
    rule = 'majority';
end
if (nargin < 6)
    validation = 'tenfold';
end

    function p = knnImplementation(x, y, l )
        [knnIdx, D] = knnsearch(x, y,'K', k,'Distance', distance);
        nnrLabels = l(knnIdx);
        if sum(nnrLabels) == 0 
            p = 0;
        elseif sum(nnrLabels) == k
            p = 1;
        elseif strcmp(rule, 'complex vote')
            p = (4 * log(min(D(~nnrLabels))) + log(mean(D(~nnrLabels))) + log(prod(D(~nnrLabels)))) >  (4 * log(min(D(nnrLabels))) + log(mean(D(nnrLabels))) + log(prod(D(nnrLabels))));
        elseif strcmp(rule, 'weighted majority')
            p = sum(1 ./ D(~nnrLabels)) < sum( 1 ./ D(nnrLabels)); 
        else
            p = round(mean(nnrLabels));
        end
    end


if strcmp(validation, 'LeaveMOut')
    r = zeros(1, 3);
    for n = 1:length(labels)
        trainIdx = true(size(labels));
        testIdx = false(size(labels));
        trainIdx(n) = false;
        testIdx(n) = true;
        testLabel = labels(testIdx);
        trainLabels = labels(trainIdx);
        train = data(trainIdx, :);
        test = data(testIdx, :);
       
        prediction = knnImplementation(train, test, trainLabels );
        r(1) = r(1) + double(prediction == testLabel);
        r(2) = r(2) + double(prediction ~= testLabel && testLabel == true); % type I error 
        r(3) = r(3) + double(prediction ~= testLabel && testLabel == false); % type II error 

    end
    r = r ./ length(labels);

elseif strcmp(validation, 'Kfold')
    indices = crossvalind('Kfold',length(labels),10);
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
            prediction = knnImplementation(train, test(j, :), trainLabels );
            rf(1) = rf(1) + double(prediction == testLabels(j));
            rf(2) = rf(2) + double(prediction ~= testLabels(j) && testLabels(j) == true); % type I error 
            rf(3) = rf(3) + double(prediction ~= testLabels(j) && testLabels(j) == false); % type II error 
        end
        r = r + rf ./ (sum(testIdx) * 10);

     end
    
else
    error('Not implemented')
end
r = r * 100;
[a, falsePositivesRate, falseNegativesRate] = deal(r(1), r(2), r(3));
% fprintf('%d-Nearest Neighbor classification results \nfor %s rule and %s distance: %.3f%% (%s)\n', k, rule, distance, a, validation);
% fprintf('False positive rate is %.3f%% and False negative rate is %.3f%%\n\n', falsePositivesRate, falseNegativesRate);

end
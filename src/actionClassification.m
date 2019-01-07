% Classify malignancy
method = 'knn';

version = 'estimated';

classificationError = struct('Accuracy', [], 'TypeI', [], 'TypeII', [], 'Input', {}, 'Projection', {}, 'Validation', {}, 'VoteRule', {}, 'Neighbours', [], 'Distance', {});
validation = {'LeaveMOut', 'Kfold'};
validationShort = {'1O', 'Kf'};
votingRules = {'majority', 'weighted majority', 'complex vote'};
votingRulesAbbr = {'M', 'W', 'C'};
distances = {'correlation', 'chebychev', 'euclidean'};
distancesAbbr = {'Cr', 'Ch', 'Eu'};
groups = {'unique', 'fixed', 'unfixed'};
projections = {'PCA', 'LDA', 'PCA->LDA', 'Spectrum'};
options.saveOptions.saveInHQ = true;


for j = 1:length(validation)
    classificationError = struct('Accuracy', [], 'TypeI', [], 'TypeII', [], 'Input', {}, 'Projection', {}, 'Validation', {}, 'VoteRule', {}, 'Neighbours', [], 'Distance', {});
    m = 0;
    for g = 1:length(groups)
        for p = 1:length(projections)
            [Gun, labels] = classifierInput(version, groups{g}, projections{p}, name);
            for i = 1:length(votingRules)
                for k = [1, 3, 5]
                    for d = 1:length(distances)
                        
                        switch method
                            case 'knn'
                                [a, b, c] = knn(Gun, labels, k, distances{d}, votingRules{i}, validation{j}); % labels = 1 for cancer

                            otherwise
                                error('Not implemented classification method');
                        end
                        
                        m = m + 1;
                        classificationError(m) = struct('Accuracy', a, 'TypeI', b, 'TypeII', c, 'Input', groups{g}, 'Projection', projections{p}, ...
                            'Validation', validation{j}, 'VoteRule', votingRules{i}, 'Neighbours', k, 'Distance', distances{d});
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
    options.saveOptions.plotName = generateName(options, ['Classification error of ', version, ' spectra with ', validation{j}]);
    plots('classificationErrors', 2, [], '', 'errors', classificationError, 'saveOptions', options.saveOptions)
end
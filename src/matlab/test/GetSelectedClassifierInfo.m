function [misclassified, falsePositives, falseNegatives] = GetSelectedClassifierInfo(ID, criterion,selectedClassifier)

    if strcmp(criterion, 'unique')
        [~, subIdx, ~] = unique(strcat({ID.SpectrumFile}, {ID.T}), 'last');

    elseif strcmp(criterion, 'unfixed')
        [~, unIdx, ~] = unique(strcat({ID.SpectrumFile}, {ID.T}), 'last');
        subIdx = intersect(unIdx, find([ID.IsFixed] == false));

    elseif strcmp(criterion, 'good')
        [~, unIdx, ~] = unique(strcat({ID.SpectrumFile}, {ID.T}), 'last');
        goodIdx = unIdx & logical([ID.IsGood]);
        %goodIdx = union(find(strcmp([ID.Sample], '9913')), union(find(strcmp([ID.Sample], '9933')), union(find(strcmp([ID.Sample], '9940')), find(strcmp([ID.Sample], '9956')))));
        subIdx = intersect(unIdx, goodIdx);

    elseif strcmp(criterion, 'all')
        subIdx = A == A;

    elseif strcmp(criterion, 'fixed')
        [~, unIdx, ~] = unique(strcat({ID.SpectrumFile}, {ID.T}), 'last');
        subIdx = intersect(unIdx, find([ID.IsFixed] == true));
        
    elseif strcmp(criterion, 'unfixedleft')
        [~, unIdx, ~] = unique(strcat({ID.SpectrumFile}, {ID.T}), 'first');
        subIdx = intersect(unIdx, find([ID.IsFixed] == false));
        
    elseif strcmp(criterion, 'unfixedright')
        [~, unIdx, ~] = unique(strcat({ID.SpectrumFile}, {ID.T}), 'last');
        subIdx = intersect(unIdx, find([ID.IsFixed] == false));
        
    else
        error('Not implemented yet.')
    end
    
    idd = ID(subIdx);

    misclassified = {idd(~selectedClassifier.Performance.IsTruePrediction).SpectrumFile}';
    falsePositives = {idd(~selectedClassifier.Performance.IsTruePrediction & contains(selectedClassifier.Performance.Labels{:}, 'Malignant')).SpectrumFile}';
    falseNegatives = {idd(~selectedClassifier.Performance.IsTruePrediction & contains(selectedClassifier.Performance.Labels{:}, 'Benign')).SpectrumFile}';
    
end
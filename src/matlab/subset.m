function [Gs, lineNamesGs, subIdx, As] = subset(input, name, criterion)

    if strcmp(input, 'measured')
        load(fullfile('..', '..', 'input', name, 'in.mat'), 'Spectra');
        G = Spectra; % G rows are observations and columns are variables

    elseif strcmp(input, 'estimated')
        load(fullfile('..', '..', 'output', name, 'ReflectanceEstimationPreset', 'out.mat'), 'EstimatedSpectra');
        G = EstimatedSpectra;

    else
        error('Not acceptable input. Choose "measured" or "estimated".')
    end

    load(fullfile('..', '..', 'input', name, 'ID.mat'), 'ID'); 

    A = [ID.IsNormal];
    X = {'Malignant', 'Benign'};
    B = X(1+A);

    lineNames = strcat([ID.Sample], '_', [ID.Type], '_', B)';

    if strcmp(criterion, 'unique')
        [~, subIdx, ~] = unique(strcat({ID.Csvid}, {ID.T}), 'last');

    elseif strcmp(criterion, 'unfixed')
        [~, unIdx, ~] = unique(strcat({ID.Csvid}, {ID.T}), 'last');
        subIdx = intersect(unIdx, find([ID.IsFixed] == false));

    elseif strcmp(criterion, 'goodright')
        [~, unIdx, ~] = unique(strcat({ID.Csvid}, {ID.T}), 'last');
        goodIdx = intersect(unIdx , find([ID.IsGood]));
        %goodIdx = union(find(strcmp([ID.Sample], '9913')), union(find(strcmp([ID.Sample], '9933')), union(find(strcmp([ID.Sample], '9940')), find(strcmp([ID.Sample], '9956')))));
        subIdx = intersect(unIdx, goodIdx);
    
    elseif strcmp(criterion, 'goodleft')
        [~, unIdx, ~] = unique(strcat({ID.Csvid}, {ID.T}), 'first');
        goodIdx = intersect(unIdx , find([ID.IsGood]));
        %goodIdx = union(find(strcmp([ID.Sample], '9913')), union(find(strcmp([ID.Sample], '9933')), union(find(strcmp([ID.Sample], '9940')), find(strcmp([ID.Sample], '9956')))));
        subIdx = intersect(unIdx, goodIdx);
        
    elseif strcmp(criterion, 'all')
        subIdx = A == A;

    elseif strcmp(criterion, 'fixed')
        [~, unIdx, ~] = unique(strcat({ID.Csvid}, {ID.T}), 'last');
        subIdx = intersect(unIdx, find([ID.IsFixed] == true));
        
    elseif strcmp(criterion, 'unfixedleft')
        [~, unIdx, ~] = unique(strcat({ID.Csvid}, {ID.T}), 'first');
        subIdx = intersect(unIdx, find([ID.IsFixed] == false));
        
    elseif strcmp(criterion, 'unfixedright')
        [~, unIdx, ~] = unique(strcat({ID.Csvid}, {ID.T}), 'last');
        subIdx = intersect(unIdx, find([ID.IsFixed] == false));
        
    else
        error('Not implemented yet.')
    end

    Gs = G(subIdx, :);
    As = A(subIdx);
    lineNamesGs = lineNames(subIdx);

end

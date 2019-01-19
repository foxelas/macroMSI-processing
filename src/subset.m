function [Gs, lineNamesGs, subIdx, As] = subset(input, name, criterion)

    if strcmp(input, 'measured')
        e = matfile(fullfile('..', '..', 'input', name, 'in.mat'));
        wavelengthN = 401;
        spectrumStruct = e.MeasuredSpectrumStruct;

    elseif strcmp(input, 'estimated')
        e = matfile(fullfile('..', '..', 'output', name, 'reflectanceestimationpreset', 'out.mat'));
        wavelengthN = 81;
        spectrumStruct = e.EstimatedSpectrumStruct;

    else
        error('Not acceptable input. Choose "measured" or "estimated".')
    end

    load(fullfile('..', '..', 'input', name, 'ID.mat'), 'ID');

    msiN = length(spectrumStruct);
    G = zeros(msiN, wavelengthN); % G rows are observations and columns are variables
    for k = 1:msiN
        G(k, :) = spectrumStruct(k).Spectrum;
    end
    A = [ID.IsNormal];
    X = {'Malignant', 'Benign'};
    B = X(1+A);

    lineNames = strcat([ID.Sample], '_', [ID.Type], '_', B)';

    if strcmp(criterion, 'unique')
        [~, subIdx, ~] = unique(strcat({ID.Csvid}, {ID.T}), 'last');

    elseif strcmp(criterion, 'unfixed')
        [~, unIdx, ~] = unique(strcat({ID.Csvid}, {ID.T}), 'last');
        subIdx = intersect(unIdx, find([ID.IsFixed] == false));

    elseif strcmp(criterion, 'good')
        [~, unIdx, ~] = unique(strcat({ID.Csvid}, {ID.T}), 'last');
        goodIdx = union(find(strcmp([ID.Sample], '9913')), union(find(strcmp([ID.Sample], '9933')), union(find(strcmp([ID.Sample], '9940')), find(strcmp([ID.Sample], '9956')))));
        subIdx = intersect(unIdx, goodIdx);

    elseif strcmp(criterion, 'all')
        subIdx = A == A;

    elseif strcmp(criterion, 'fixed')
        [~, unIdx, ~] = unique(strcat({ID.Csvid}, {ID.T}), 'last');
        subIdx = intersect(unIdx, find([ID.IsFixed] == true));

    else
        error('Not implemented yet.')
    end

    Gs = G(subIdx, :);
    As = A(subIdx);
    lineNamesGs = lineNames(subIdx);

end

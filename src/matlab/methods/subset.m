function [inputSubset, subsetSamples, subsetIndexes, malignancy] = subset(input, name, criterion)
%%SUBSET select a subset of the dataset based on a criterion 
% 
% Inputs: 
% input - dataset for input (measured or estimated)
% name - dataset version name 
% criterion - 'unique', 'unfixed', 'fixed', 'all', etc.
% 
% Outputs: 
% inputSubset - the subset
% subsetSamples - the names of the subset samples 
% subsetIndexes - indexes of the samples contained in the subset 
% malignancy - true/false values for malignancy in the subset
% 
% Usage: 
% [inputSubset, subsetSamples, subsetIndexes, malignancy] = subset('estimated', 'saitama_v3', 'fixed')
% 
    if strcmp(input, 'measured')
        load(fullfile(generateName('input'), name, 'in.mat'), 'Spectra');
        G = Spectra; % G rows are observations and columns are variables

    elseif strcmp(input, 'estimated')
        load( fullfile(generateName( 'output'), name, 'ReflectanceEstimationPreset', 'out.mat'), 'EstimatedSpectra');
        G = EstimatedSpectra;

    else
        error('Not acceptable input. Choose "measured" or "estimated".')
    end

    load(fullfile(generateName('input'), name, 'ID.mat'), 'ID'); 

    benignity = [ID.IsBenign];
    X = {'Malignant', 'Benign'};
    B = X(1+benignity);

    sampleNames = strcat([ID.Sample], '_', [ID.Type], '_', B)';

    if strcmp(criterion, 'unique')
        [~, subsetIndexes, ~] = unique(strcat({ID.SpectrumFile}, {ID.T}), 'last');

    elseif strcmp(criterion, 'unfixed')
        [~, unIdx, ~] = unique(strcat({ID.SpectrumFile}, {ID.T}), 'last');
        subsetIndexes = intersect(unIdx, find([ID.IsFixed] == false));

    elseif strcmp(criterion, 'goodright')
        [~, unIdx, ~] = unique(strcat({ID.SpectrumFile}, {ID.T}), 'last');
        goodIdx = intersect(unIdx , find([ID.IsGood]));
        %goodIdx = union(find(strcmp([ID.Sample], '9913')), union(find(strcmp([ID.Sample], '9933')), union(find(strcmp([ID.Sample], '9940')), find(strcmp([ID.Sample], '9956')))));
        subsetIndexes = intersect(unIdx, goodIdx);
    
    elseif strcmp(criterion, 'goodleft')
        [~, unIdx, ~] = unique(strcat({ID.SpectrumFile}, {ID.T}), 'first');
        goodIdx = intersect(unIdx , find([ID.IsGood]));
        %goodIdx = union(find(strcmp([ID.Sample], '9913')), union(find(strcmp([ID.Sample], '9933')), union(find(strcmp([ID.Sample], '9940')), find(strcmp([ID.Sample], '9956')))));
        subsetIndexes = intersect(unIdx, goodIdx);
        
    elseif strcmp(criterion, 'all')
        subsetIndexes = benignity == benignity;

    elseif strcmp(criterion, 'fixed')
        [~, unIdx, ~] = unique(strcat({ID.SpectrumFile}, {ID.T}), 'last');
        subsetIndexes = intersect(unIdx, find([ID.IsFixed] == true));
        
    elseif strcmp(criterion, 'unfixedleft')
        [~, unIdx, ~] = unique(strcat({ID.SpectrumFile}, {ID.T}), 'first');
        subsetIndexes = intersect(unIdx, find([ID.IsFixed] == false));
        
    elseif strcmp(criterion, 'unfixedright')
        [~, unIdx, ~] = unique(strcat({ID.SpectrumFile}, {ID.T}), 'last');
        subsetIndexes = intersect(unIdx, find([ID.IsFixed] == false));
        
    else
        error('Not implemented yet.')
    end

    inputSubset = G(subsetIndexes, :);
    malignancy = ~benignity(subsetIndexes);
    subsetSamples = sampleNames(subsetIndexes);

end

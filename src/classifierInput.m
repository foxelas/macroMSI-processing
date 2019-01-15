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

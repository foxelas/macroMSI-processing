function [G, labels] = classifierInput(version, group, projection, name, labelsAsText)

if (nargin < 5)
    labelsAsText = false;
end

[Gfx, ~, ~, labelsfx] = subset(version, name, group);

if contains(projection, 'pca' ) || contains(projection, 'lda')
    [~, scores] = dimensionReduction(projection, Gfx, double(labelsfx));
    G = scores(:, 1:10);
    
elseif strcmp(projection, 'spectrum')
    G = Gfx;
    
else 
    G = Gfx;
end

labels = ~labelsfx';

if (labelsAsText)
    X = {'Benign', 'Malignant'};
    labels = X(1+labels);
end

end

color = outfile.EstimatedSpectra';
mmlbp = outfile.MMLbpFeatures;
texture = [mmlbp{1}, mmlbp{2}]';

input = color';
% input = [color; texture]';
labels = {ID.Label};
%labels = oldLabels;

dimred = 'LDA';
[WMeasured, score, latent, explained] = reduceDimension(dimred, input, labels);
plotDimensionReduction(strcat(dimred, 'color'), strcat('2 Class ', dimred), score, labels, latent(1:10), explained(1:10), 1, options.saveOptions);

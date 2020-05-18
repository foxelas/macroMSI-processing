outfile = matfile(fullfile(options.saveOptions.savedir, outputFolderMap('features'), 'out.mat'));
color = outfile.EstimatedSpectra';
mmlbp = outfile.MMLbpFeatures;
texture = [mmlbp{1}, mmlbp{2}]';

%% only color
inputs = color;
dimension1 = 5;
dimension2 = 5;
[~, outputs] = getSom(inputs, dimension1, dimension2, options.saveOptions, 'spect');
plotAssignToNeurons(outputs, {ID.Label});

%% color and texture
inputs = [color; texture];

% Create a Self-Organizing Map
dimension1 = 10;
dimension2 = 10;
[~, outputs] = getSom(inputs, dimension1, dimension2, options.saveOptions, 'spect+mmlbp');
plotAssignToNeurons(outputs, {ID.Label});


function [] = plotAssignToNeurons(outputs, labels)
clusters = zeros(size(outputs));
for i = 1:size(outputs, 2)
    if strcmp(labels{i}, 'Atypical')
        clusters(:, i) = outputs(:, i) * 2;
    elseif strcmp(labels{i}, 'Malignant')
        clusters(:, i) = outputs(:, i) * 3;
    else
        clusters(:, i) = outputs(:, i);
    end
end
figure(5);
imagesc(clusters);
colormap('jet');
colorbar;
end
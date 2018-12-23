function [G, labels] = classifierInput(version, group, projection, name)
    
[Gfx, ~, ~, labelsfx] = subset(version, name, group);

if strcmp(projection, 'PCA')
    [~, scores, ~, ~, ~]= pca(Gfx,'Centered', true);   
    G = scores(:, 1:10);
elseif strcmp(projection, 'LDA')
    [~, scores] = LDA(Gfx, double(labelsfx), [], 'pdf');
    G = scores(:, 1:10);

elseif strcmp(projection, 'PCALDA')
    [~, scores, ~, ~, ~]= pca(Gfx,'Centered', true);   
    [~, G] = LDA(scores(:, 1:10), double(labelsfx), [], 'pdf');

else
    G = Gfx;
end

labels = ~labelsfx';

end


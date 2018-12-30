function [W, score, latent, explained] = dimensionReduction(method, input, labels, additionalMethod, priors)
    latent = [];
    explained = [];
    
    switch method
        case 'pca'
            [W, score, latent, explained] = PCA(input, 'Centered', true);
            
        case 'lda'
            [W, score] = LDA(input, labels, [], additionalMethod);
            
        case 'pcalda'
            [~, score, ~, ~] = PCA(input, 'Centered', true);
            [W, score] = LDA(score(:, 1:10), labels, [], additionalMethod);
            
        otherwise
            error('Not implemented dimension reduction method')
    end
end


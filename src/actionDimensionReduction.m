%Data projection for dimension reduction (pca, lda, etc)
additionalMethod = ''; 

if contains(action, 'pcalda')
    dimredMethod = 'pcalda';
    
elseif contains(action, 'pca b')
    dimredMethod = 'pca b';
    
elseif contains(action, 'lda b')
    dimredMethod = 'lda b';
    
elseif contains(action, 'lda')
    dimredMethod = 'lda';
    
elseif contains(action, 'pca')
    dimredMethod = 'pca';
    
else 
    error('Unsupported dimension reduction method');
end

input = 'estimated';
[Gun, lineNamesun, ~, labelsun] = subset(input, name, 'unique');
[W, score, latent, explained] = dimensionReduction(dimredMethod, Gun, double(labelsun));
%L = [ones(length(un),1) Gun] * W';
%P = exp(L) ./ repmat(sum(exp(L),2),[1 2]);
options.saveOptions.plotName = generateName(options, [dimredMethod, ' of ', input, ' spectra by type']);
plots(dimredMethod, 1, score, [dimredMethod,' Fix'], 'LineNames', lineNamesun, 'Latent', latent(1:10), 'Explained', explained(1:10), 'SaveOptions', options.saveOptions)
options.saveOptions.plotName = generateName(options, [dimredMethod, ' of ', input, ' spectra by sample']);
plots(dimredMethod, 2, score, [dimredMethod, ' Sample'], 'LineNames', lineNamesun, 'Latent', latent(1:10), 'Explained', explained(1:10), 'SaveOptions', options.saveOptions)

%unfixed data only
[Gfx, lineNamesfx, ~, labelsfx] = subset(input, name, 'unfixed');
[W, score, latent, explained] = dimensionReduction(dimredMethod, Gfx, double(labelsfx));
options.saveOptions.plotName = generateName(options, [dimredMethod, ' of ', input, ' spectra by sample (only unfixed)']);
plots(dimredMethod, 3, score, [dimredMethod,' Sample'], 'LineNames', lineNamesfx, 'Latent', latent(1:10), 'Explained', explained(1:10),'SaveOptions', options.saveOptions)

input = 'measured';
[Gun, lineNamesun, ~, labelsun] = subset(input, name, 'unique');
[W, score, latent, explained] = dimensionReduction(dimredMethod, Gun, double(labelsun));
%L = [ones(length(un),1) Gun] * W';
%P = exp(L) ./ repmat(sum(exp(L),2),[1 2]);
options.saveOptions.plotName = generateName(options, [dimredMethod, ' of ', input, ' spectra by type']);
plots(dimredMethod, 4, score, [dimredMethod,' Fix'], 'LineNames', lineNamesun, 'Latent', latent(1:10), 'Explained', explained(1:10),'SaveOptions', options.saveOptions)
options.saveOptions.plotName = generateName(options, [dimredMethod, ' of ', input, ' spectra by sample']);
plots(dimredMethod, 5, score, [dimredMethod, ' Sample'], 'LineNames', lineNamesun, 'Latent', latent(1:10), 'Explained', explained(1:10),'SaveOptions', options.saveOptions)

%unfixed data only
[Gfx, lineNamesfx, ~, labelsfx] = subset(input, name, 'unfixed');
[W, score, latent, explained] = dimensionReduction(dimredMethod, Gfx, double(labelsfx));
options.saveOptions.plotName = generateName(options, [dimredMethod,' of ', input, ' spectra by sample (only unfixed)']);
plots(dimredMethod, 6, score, [dimredMethod,' Sample'], 'LineNames', lineNamesfx, 'Latent', latent(1:10), 'Explained', explained(1:10),'SaveOptions', options.saveOptions)


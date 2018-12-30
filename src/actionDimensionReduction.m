additionalMethod = ''; 

if contains(action, 'pca')
    dimredMethod = 'pca';
elseif contains(action, 'lda')
    dimredMethod = 'lda';
    additionalMethod = 'p';  % 'pdf', 'kardi', ''
elseif contains(action, 'pcalda')
    dimredMethod = 'pcalda';
    additionalMethod = 'p';  % 'pdf', 'kardi', ''
else 
    dimredMethod = 'pca';
end

input = 'estimated';
[Gun, lineNamesun, ~, labelsun] = subset(input, name, 'unique');
[W, score, latent, explained] = dimensionReduction(dimredMethod, Gun, double(labelsun), additionalMethod, []);
%L = [ones(length(un),1) Gun] * W';
%P = exp(L) ./ repmat(sum(exp(L),2),[1 2]);
options.saveOptions.plotName = generateName(options, [dimredMethod, ' of ', input, ' spectra by type', ldaMethod]);
plots(dimredMethod, 1, score, [dimredMethod,' Fix'], 'lineNames', lineNamesun, 'saveOptions', options.saveOptions)
options.saveOptions.plotName = generateName(options, [dimredMethod, ' of ', input, ' spectra by sample', ldaMethod]);
plots(dimredMethod, 2, score, [dimredMethod, ' Sample'], 'lineNames', lineNamesun, 'latent', latent(1:10), 'explained', explained(1:10), 'saveOptions', options.saveOptions)

%unfixed data only
[Gfx, lineNamesfx, ~, labelsfx] = subset('estimated', name, 'unfixed');
[W, score, latent, explained] = dimensionReduction(dimredMethod, Gfx, double(labelsfx),  additionalMethod, []);
options.saveOptions.plotName = generateName(options, [dimredMethod, ' of ', input, ' spectra by sample (only unfixed)', ldaMethod]);
plots(dimredMethod, 3, score, [dimredMethod,' Sample'], 'lineNames', lineNamesfx, 'latent', latent(1:10), 'explained', explained(1:10),'saveOptions', options.saveOptions)

input = 'measured';
[Gun, lineNamesun, ~, labelsun] = subset(input, name, 'unique');
[W, score, latent, explained] = dimensionReduction(dimredMethod, Gun, double(labelsun),  additionalMethod, []);
%L = [ones(length(un),1) Gun] * W';
%P = exp(L) ./ repmat(sum(exp(L),2),[1 2]);
options.saveOptions.plotName = generateName(options, [dimredMethod, ' of ', input, ' spectra by type', ldaMethod]);
plots(dimredMethod, 4, score, [dimredMethod,' Fix'], 'lineNames', lineNamesun, 'latent', latent(1:10), 'explained', explained(1:10),'saveOptions', options.saveOptions)
options.saveOptions.plotName = generateName(options, [dimredMethod, ' of ', input, ' spectra by sample', ldaMethod]);
plots(dimredMethod, 5, score, [dimredMethod, ' Sample'], 'lineNames', lineNamesun, 'latent', latent(1:10), 'explained', explained(1:10),'saveOptions', options.saveOptions)

%unfixed data only
[Gfx, lineNamesfx, ~, labelsfx] = subset('estimated', 'unfixed', options.saveOptions.savedir);
[W, score, latent, explained] = dimensionReduction(dimredMethod, Gfx, double(labelsfx),  additionalMethod, []);
options.saveOptions.plotName = generateName(options, [dimredMethod,' of ', input, ' spectra by sample (only unfixed)', ldaMethod]);
plots(dimredMethod, 6, score, [dimredMethod,' Sample'], 'lineNames', lineNamesfx, 'latent', latent(1:10), 'explained', explained(1:10),'saveOptions', options.saveOptions)


%Data projection for dimension reduction (PCA, LDA, etc)
additionalMethod = ''; 
action = options.action; 

if contains(action, 'PCALDA')
    dimredMethod = 'PCALDA';
    
elseif contains(action, 'PCA b')
    dimredMethod = 'PCA b';
    
elseif contains(action, 'LDA b')
    dimredMethod = 'LDA b';
    
elseif contains(action, 'LDA')
    dimredMethod = 'LDA';
    
elseif contains(action, 'PCA')
    dimredMethod = 'PCA';
    
else 
    error('Unsupported dimension reduction method');
end
dimredMethodCap = upper(dimredMethod);

w = warning('on', 'all');

%% Measured
input = 'measured';
[Gun, lineNamesun, ~, labelsun] = subset(input, name, 'unique');
[WMeasured, score, latent, explained] = dimensionReduction(dimredMethod, Gun, double(labelsun));
%L = [ones(length(un),1) Gun] * W';
%P = exp(L) ./ repmat(sum(exp(L),2),[1 2]);
options.saveOptions.plotName = generateName([dimredMethodCap, ' of ', input, ' spectra by type'], options);
plots(dimredMethod, 1, score, [dimredMethodCap,' Fix'], 'LineNames', lineNamesun, 'Latent', latent(1:10), 'Explained', explained(1:10),'SaveOptions', options.saveOptions)
options.saveOptions.plotName = generateName([dimredMethodCap, ' of ', input, ' spectra by sample'], options);
plots(dimredMethod, 2, score, [dimredMethodCap, ' Sample'], 'LineNames', lineNamesun, 'Latent', latent(1:10), 'Explained', explained(1:10),'SaveOptions', options.saveOptions)

%unfixed data only
[Gfx, lineNamesfx, ~, labelsfx] = subset(input, name, 'unfixed');
[WMeasuredUnfixed, score, latent, explained] = dimensionReduction(dimredMethod, Gfx, double(labelsfx));
options.saveOptions.plotName = generateName([dimredMethodCap,' of ', input, ' spectra by sample (only unfixed)'], options);
plots(dimredMethod, 3, score, [dimredMethodCap,' Sample'], 'LineNames', lineNamesfx, 'Latent', latent(1:10), 'Explained', explained(1:10),'SaveOptions', options.saveOptions)

%% Estimated spectra projected by measured projection base 
options.saveOptions.plotName = generateName([dimredMethodCap , ' of ', input, ' spectra by type [by projection base of measured]'], options);
plots(dimredMethod, 7, Gun * WMeasured, [dimredMethodCap,' Fix'], 'LineNames', lineNamesun, 'SaveOptions', options.saveOptions);
options.saveOptions.plotName = generateName([dimredMethodCap, ' of ', input, ' spectra by sample [by projection base of measured]'], options);
plots(dimredMethod, 8, Gun * WMeasured, [dimredMethodCap, ' Sample'], 'LineNames', lineNamesun, 'SaveOptions', options.saveOptions)

options.saveOptions.plotName = generateName([dimredMethodCap, ' of ', input, ' spectra by sample (only unfixed) [by projection base of measured]'], options);
plots(dimredMethod, 9, Gfx * WMeasuredUnfixed, [dimredMethodCap,' Sample'], 'LineNames', lineNamesfx, 'SaveOptions', options.saveOptions)

warning(w)
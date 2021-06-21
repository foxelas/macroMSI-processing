%Data projection for dimension reduction (PCA, LDA, etc)
additionalMethod = '';
action = getSetting('action');

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
[WMeasured, score, latent, explained] = reduceDimension(dimredMethod, Gun, double(labelsun));
%L = [ones(length(un),1) Gun] * W';
%P = exp(L) ./ repmat(sum(exp(L),2),[1 2]);
setSetting('plotName', generateName([dimredMethodCap, ' of ', input, ' spectra by type']));
%plots(dimredMethod, 1, score, [dimredMethodCap, ' Fix'], 'LineNames', lineNamesun, 'Latent', latent(1:10), 'Explained', explained(1:10))
setSetting('plotName', generateName([dimredMethodCap, ' of ', input, ' spectra by sample']));
%plots(dimredMethod, 2, score, [dimredMethodCap, ' Sample'], 'LineNames', lineNamesun, 'Latent', latent(1:10), 'Explained', explained(1:10))

%unfixed data only
[Gfx, lineNamesfx, ~, labelsfx] = subset(input, name, 'unfixed');
[WMeasuredUnfixed, score, latent, explained] = reduceDimension(dimredMethod, Gfx, double(labelsfx));
setSetting('plotName', generateName([dimredMethodCap, ' of ', input, ' spectra by sample (only unfixed)']));
%plots(dimredMethod, 3, score, [dimredMethodCap, ' Sample'], 'LineNames', lineNamesfx, 'Latent', latent(1:10), 'Explained', explained(1:10))

%% Estimated spectra projected by measured projection base
setSetting('plotName', generateName([dimredMethodCap, ' of ', input, ' spectra by type [by projection base of measured]']));
%plots(dimredMethod, 7, Gun*WMeasured, [dimredMethodCap, ' Fix'], 'LineNames', lineNamesun);
setSetting('plotName', generateName([dimredMethodCap, ' of ', input, ' spectra by sample [by projection base of measured]']));
%plots(dimredMethod, 8, Gun*WMeasured, [dimredMethodCap, ' Sample'], 'LineNames', lineNamesun)

setSetting('plotName', generateName([dimredMethodCap, ' of ', input, ' spectra by sample (only unfixed) [by projection base of measured]']));
plots(fig, @plotDimensionReduction, dimredMethod, [dimredMethodCap, ' Sample'], Gfx*WMeasuredUnfixed);
%plots(dimredMethod, 9, Gfx*WMeasuredUnfixed, [dimredMethodCap, ' Sample'], 'LineNames', lineNamesfx);

warning(w)
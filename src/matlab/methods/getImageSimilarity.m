function [metrics] = getImageSimilarity(base, ref, filename, id1, id2, name1, name2)
%   Usage
%   [metrics] = getImageSimilarity(pc_unfixed, pc_fixed, strcat('_pc', ...
%   num2str(i)), getSetting('pca'), unfixedId, fixedId, name1, name2);

savedir = getSetting('savedir');
if ~isempty(filename)
    filename = fullfile(savedir, '14-MapComparison', filename, strcat(num2str(id1), 'vs', num2str(id2), '_', filename));
end

% setSetting('plotName', strcat(filename, '.png'));
% plotFunWrapper(1, @plotMontage, base, ref, strcat(name1, ' vs ', name2));

baseMask = ~isnan(base);
refMask = ~isnan(ref);
minVal = 0.99 * min(min(base(:)), min(ref(:)));
baseWithoutNan = base;
baseWithoutNan(~baseMask) = minVal;
refWithoutNan = ref;
refWithoutNan(~refMask) = minVal;

%% SSIM
[ssimval, ssimmap] = ssim(baseWithoutNan, refWithoutNan);
fprintf('Global SSIM Value is %.5f \n', ssimval);
%{
figTitle = strcat('Local SSIM Map with Global SSIM Value: ',num2str(ssimval));
setSetting('plotName', strcat(filename, '_ssim.png'));
plotFunWrapper(2, @plotScaledImage, ssimmap, figTitle, []);
%}

%% Sum of squared differences (SSD)
ssd = sum((baseWithoutNan(:) - refWithoutNan(:)).^2);
fprintf('Sum of squared differences is %.5f \n', ssd);

%% Normalized cross correlation coefficient (NCC)
covariance = cov(base, ref, 'omitrows');
ncc = covariance(1, 2)^2 / (var(base(:), 'omitnan') * var(ref(:), 'omitnan')); %'partialrows'
fprintf('Normalized cross correlation coefficient is %.5f \n', ncc);

%% Pearsons correlation coefficient
cc = corrcoef(base, ref, 'Rows', 'pairwise');
cc = cc(1, 2);
fprintf('Pearsons Correlation Coefficients is %.5f \n', cc);
%{
ccmap = xcorr2(base,ref);
figTitle = strcat('Pearsons Correlation Coefficients: ',num2str(ccmap));
setSetting('plotName', strcat(filename, '_cc.png'));
plotFunWrapper(3, @plotScaledImage, ccmap, figTitle, []);
%}

%% Correlation Coefficient Ration
cr = nan;

%% Similarity of histogram intersection
n = 20;
figure(3);
hBase = histogram(base(baseMask), n);
setSetting('plotName', strcat(filename, '_hist_', name1, '.png'));
savePlot(gcf);
counts1 = hBase.Values;
edges = hBase.BinEdges;
centers1 = edges(1:end-1)' + diff(edges)' / 2;
figure(4);
hRef = histogram(ref(refMask), n);
setSetting('plotName', strcat(filename, '_hist_', name2, '.png'));
savePlot(gcf);
counts2 = hRef.Values;
edges = hRef.BinEdges;
centers2 = edges(1:end-1)' + diff(edges)' / 2;
similarityHistIntersection = sum(min(counts1, counts2)) / sum(counts2);
fprintf('Similarity of histogram intersection is %.5f \n', similarityHistIntersection);

%% Kullback Leiler divergence
klDist = KLDiv(counts2, counts1); %P||Q from Q to P , Q is prior, P is posterior
fprintf('Kullback-Leibler divergence from Unfixed to Fixed of histograms is %.5f \n', klDist);

%% Earth Mover's Distance
w1 = counts1' / sum(counts1);
w2 = counts2' / sum(counts2);
[flow, fval] = emd(centers1, centers2, w1, w2, @gdf);
fprintf('Earth Movers Distance of histograms is %.5f \n', fval);

%% Results
metrics = [ssimval, ssd, ncc, cc, cr, similarityHistIntersection, klDist, fval];
% metricsMaps = {ssimmap, ccmap};

end

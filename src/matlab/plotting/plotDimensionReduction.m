function [] = plotDimensionReduction(dimred, figTitle, coefficients, labels, latent, explained, fig,saveOptions)

    if (nargin < 5)
		latent = [];
    end
    if (nargin < 6)
		explained = [];
    end
    if (nargin < 7)
        fig = figure;
    else     
        figure(fig);
        clf(fig);
    end
    if (nargin < 8)
        saveOptions.SaveImage = false;
    end
  
	
	if strcmp(dimred, 'PCA') && ~isempty(latent) && ~isempty(explained)
		subplot(1, 3, 1);
		plot(1:length(latent), latent, '-mx');
		xlabel('Sorted eigenvalue index');
		ylabel('PCA eigenvalues');
		xlim([1, length(latent) + 1])
		text(1:length(explained), latent+10^(ceil(log10(max(latent))) - 2), arrayfun(@(x) sprintf('%.2f%%', x), explained, 'UniformOutput', false), 'FontSize', 8);
		title('Largest eigenvalues and explained covariance percentage')
		
		subplot(1, 3, [2, 3]);
	end
	
	marker = ['o', 'x', 'd', '^', '*', 'h', 'p', 'v', 's', '<', '+', '>'];
	observations = size(coefficients, 1);
	attr = split(labels, ' ');
	sample = {attr{:, 1}}';
	type = {attr{:, 2}}';
	isBenign = {attr{:, 3}}';
	idx = strcmp(isBenign, 'Normal') | strcmp(isBenign, 'Benign');
	colors(idx) = 'b';
	colors(~ismember(idx, 1:length(isBenign))) = 'r';
	
	%dummy legends
	h = zeros(2, 1);
	hold on
	h(1) = plot([NaN, NaN], 'Color', 'r', 'DisplayName', 'Malignant');
	h(2) = plot([NaN, NaN], 'Color', 'b', 'DisplayName', 'Benign');
	hold off
	
	if contains(figTitle, 'Fix', 'IgnoreCase', true)
		hold on
		h(3) = plot([NaN, NaN], 'Color', 'k', 'Marker', marker(1), 'DisplayName', 'Fixed');
		h(4) = plot([NaN, NaN], 'Color', 'k', 'Marker', marker(2), 'DisplayName', 'Unfixed');
		hold off
		markers = repmat(marker(1), observations, 1);
		idx = arrayfun(@(x) any(strcmp(x, 'unfixed')), type);
		markers(idx) = marker(2);
	end
	
	if contains(figTitle, 'Sample', 'IgnoreCase', true)
		markers = repmat(marker(1), observations, 1);
		samples = unique(sample, 'stable');
		for i = 1:length(samples)
			hold on
			h(i+2) = plot([NaN, NaN], 'Color', 'k', 'Marker', marker(i), 'DisplayName', samples{i});
			hold off
			idx = arrayfun(@(x) any(strcmp(x, samples{i})), sample);
			markers(idx) = marker(i);
		end
	end
	%dummy legends
	
	hold on
	for i = 1:observations
		scatter(coefficients(i, 1), coefficients(i, 2), [], colors(i), markers(i));
	end
	hold off
	
	if strcmp(dimred, 'pca')
		title('PCA scores for PC1&2')
		xlabel('Principal component 1');
		ylabel('Principal component 2');
	else
		title('LDA projections for LD1&2')
		xlabel('Lidear discriminant 1');
		ylabel('Linear Discriminant 2');
	end
	ax = gca;
	ax.XRuler.Exponent = 0;
	ax.YRuler.Exponent = 0;
	titl = strsplit(saveOptions.plotName, '\');
	figTitle = strrepAll(titl{end});
	suptitle(figTitle);
	
	legend(h, 'Location', 'best');
	set(gcf, 'Position', get(0, 'Screensize'));
	
	saveOptions.saveInHQ = true;	
    savePlot(fig, saveOptions);

end


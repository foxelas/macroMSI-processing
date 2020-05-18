function [] = plotreduceDimension(dimred, figTitle, coefficients, labels, latent, explained, fig,saveOptions)

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
  
	
	if contains(dimred, 'PCA') && ~isempty(latent) && ~isempty(explained)
		subplot(1, 3, 1);
		plot(1:length(latent), latent, '-mx');
		xlabel('Sorted eigenvalue index');
		ylabel('PCA eigenvalues');
		xlim([1, length(latent) + 1])
		text(1:length(explained), latent+10^(ceil(log10(max(latent))) - 2), arrayfun(@(x) sprintf('%.2f%%', x), explained, 'UniformOutput', false), 'FontSize', 8);
		title('Largest eigenvalues and explained covariance percentage')
		
		subplot(1, 3, [2, 3]);
	end
	
	observations = size(coefficients, 1);
    
	%dummy legends
    markerStyleMap = getMarkerStyleMap(getStyleName(labels)); 
    lineColorMap = getLineColorMap(getStyleName(labels));
    
    key = keys(lineColorMap);
    hold on; 
    h = zeros(1,length(key));
    for i = 1:length(key)
        h(i) = plot(nan, nan, markerStyleMap(key{i}), 'MarkerSize', 12,  'Color', lineColorMap(key{i}), ...
                'MarkerFaceColor', lineColorMap(key{i}),  'LineWidth', 3, 'DisplayName', key{i});
    end
    hold off;
    
    hold on
	for i = 1:observations
		scatter(coefficients(i, 1), coefficients(i, 2), [], ...
            lineColorMap(labels{i}), markerStyleMap(labels{i}), ...            
            'MarkerEdgeColor', lineColorMap(labels{i}), ... 'k',...
            'MarkerFaceColor', lineColorMap(labels{i}), ...
            'LineWidth', 3);
	end
	hold off
	
	if contains(dimred, 'PCA')
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
    
    if contains(dimred, 'PCA') && ~isempty(latent) && ~isempty(explained)
        sgtitle(figTitle);
    end
	
	legend(h, 'Location', 'best');
	set(gcf, 'Position', get(0, 'Screensize'));
	
	saveOptions.saveInHQ = true;
    outputFolderMap = getOutputDirectoryMap();
    saveOptions.plotName = fullfile(saveOptions.savedir, outputFolderMap('dimred'), dimred);
    savePlot(fig, saveOptions);

end


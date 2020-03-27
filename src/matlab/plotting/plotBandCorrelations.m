function [] = plotBandCorrelations(correlations, msiType, figTitle, fig, saveOptions)
  
	if (nargin < 4)
		fig = figure;
	else     
		figure(fig);
        clf(fig);
	end
	if (nargin < 5)
		saveOptions.SaveImage = false;
    end
    
    imagesc(correlations);
    xticklabels(getBands(msiType));
    yticklabels(getBands(msiType));
    xlabel('Bands');
    ylabel('Bands');
    %title(sprintf('Band Correlation for %s', patchNames{i}));
    title(figTitle);
    colorbar();
    caxis([0 1]);
    
    saveOptions.cropBorders = true;
	savePlot(fig, saveOptions);
    
end 
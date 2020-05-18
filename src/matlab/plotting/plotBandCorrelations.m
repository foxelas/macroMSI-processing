function [] = plotBandCorrelations(correlations, msiType, figTitle, fig)
  
    imagesc(correlations);
    xticklabels(getBands(msiType));
    yticklabels(getBands(msiType));
    xlabel('Bands');
    ylabel('Bands');
    %title(sprintf('Band Correlation for %s', patchNames{i}));
    title(figTitle);
    colorbar();
    caxis([0 1]);
    
    setSetting('cropBorders', true);
	savePlot(fig);
    
end 
function [] = plotFeatureImportance(importance,fig)
%% Plot important features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (nargin < 1)
        importance = [0.02170081 0.03423029 0.04734384 0.07019354 0.08864051 0.1096366 0.13921365 0.14391144 0.17794282];
    end
    
    features = categorical({'760', '635', '620', '650', '645', '655', '665', '625', '660'});
    bar(features, importance);
    title('Top 10 Important Wavelengths', 'FontSize', 15);
    ylabel('Relative Importance', 'FontSize', 12);
    xlabel('Wavelength (nm)', 'FontSize', 12);
    setSetting('plotName', fullfile(getSetting('savedir'),'Importance', 'featimp'));
    savePlot(fig);

end


function [] = plotMap(I, fig, saveOptions)

    if (nargin < 3)
        fig = figure;
    else 
        figure(fig);
        clf(fig);
    end
    
    c = imagesc(I);
    c.Parent.Visible = 'off';
    colorbar;
    title(saveOptions.figTitle);
    
    saveOptions.cropBorders = true;
    saveOptions.plotName = fullfile(saveOptions.savedir, getOutputDirectoryMap('map'), saveOptions.relativeDir, saveOptions.outName);

    savePlot(fig, saveOptions);

end


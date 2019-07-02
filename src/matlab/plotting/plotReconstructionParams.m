function [] = plotReconstructionParams(x, y, xlab, ylab, figTitle, fig,saveOptions)

    if (nargin < 6)
        fig = figure;
    else 
        figure(fig);
        clf(fig);
    end

    if (nargin < 7)
        saveOptions.SaveImage = false;
    end 

    warning('off');
    
    plot(log10(x), y, 'LineWidth', 5);
    xlabel(xlab, 'FontSize', 15);
    ylabel(ylab, 'FontSize', 15);
    title(figTitle, 'FontSize', 15);

    savePlot(fig, saveOptions);
    
    warning('on');
end
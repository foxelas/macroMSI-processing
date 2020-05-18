function [] = plotReconstructionParams(x, y, xlab, ylab, figTitle, fig)

    warning('off');
    
    plot(log10(x), y, 'LineWidth', 5);
    xlabel(xlab, 'FontSize', 15);
    ylabel(ylab, 'FontSize', 15);
    title(figTitle, 'FontSize', 15);

    savePlot(fig);
    
    warning('on');
end
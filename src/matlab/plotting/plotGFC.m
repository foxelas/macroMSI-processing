function [] = plotGFC(x, y, xlab, ylab, figTitle, fig,saveOptions)

    if (nargin < 7)
        fig = figure;
    else 
        figure(fig);
        clf(fig);
    end

    if (nargin < 8)
        saveOptions.SaveImage = false;
    end 

    warning('off');
    
    if contains(xlab, 'sigma')
        x = log10(x);
    end
    plot(x, y, 'LineWidth', 5);
    xlabel(xlab, 'FontSize', 15, 'Interpreter','latex');
    ylabel(ylab, 'FontSize', 15);
    %title(figTitle,'FontSize', 15);
    ylim([0.95, 1.005]);
    yline(0.99, '--r', 'Good', 'FontSize', 15, 'LabelHorizontalAlignment', 'left');
    yline(0.999, '--m', 'Very Good',  'FontSize', 15, 'LabelHorizontalAlignment', 'left');
    savePlot(fig, saveOptions);
    
    warning('on');
end
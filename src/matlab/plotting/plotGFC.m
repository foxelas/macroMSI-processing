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
    
    if (size(y,1) ~= size(x,1))
        error('Incompatible size.');
    end
    
    if contains(xlab, 'sigma')
        x = log10(x);
    end
    plot(x, y, 'LineWidth', 5);
    xlabel(xlab, 'FontSize', 15, 'Interpreter','latex');
    ylabel(ylab, 'FontSize', 15);
    %title(figTitle,'FontSize', 15);
    ylim([0.95, 1.005]);
    yline(0.99, '--b', 'Good');
    yline(0.999, '--m', 'Very Good');
    savePlot(fig, saveOptions);
    
    warning('on');
end
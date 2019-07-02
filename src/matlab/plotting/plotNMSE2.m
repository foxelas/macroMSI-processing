function [] = plotNMSE2(x1, x2, y, x1lab, x2lab, ylab, figTitle, fig,saveOptions)

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
    
    y = reshape(y, [length(x1), length(x2)]);
    xy = meshgrid(log10(x1), log10(x2));
    s = surf(xy, y, 'FaceAlpha',0.5);
    s.EdgeColor = 'none';
    xlabel(x1lab, 'FontSize', 15, 'Interpreter','latex');
    ylabel(x2lab, 'FontSize', 15, 'Interpreter','latex');
    zlabel(ylab, 'FontSize', 15);
    %title(figTitle, 'FontSize', 15);
    colorbar

    savePlot(fig, saveOptions);
    
    warning('on');
end
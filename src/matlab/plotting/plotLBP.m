function [] = plotLBP( lbpImage, figTitle, fig,saveOptions)
%% Plot LBP values on image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (nargin < 3)
        fig = figure;
    else     
        figure(fig);
        clf(fig);
    end
    if (nargin < 4)
        saveOptions.SaveImage = false;
    end
    
    im = imagesc(imrotate(lbpImage,90));
    colormap('parula');
    if max(lbpImage(:)) > 0.3
        upperlimit = 0.4;
    else 
        upperlimit = 0.04;
    end
    im.AlphaData = .8;
    c = colorbar;
    axis off;
    c.Label.FontSize = 15;
    c.Label.FontWeight = 'bold';
    c.Label.String = 'Texture Descriptor Value';
    c.Limits = [0,upperlimit];
    c.LimitsMode = 'manual';

%     title(strjoin({figTitle, 'descriptor values'}, ' '),'FontSize', 20);
%     set(gcf, 'Position', get(0, 'Screensize'));
    
    %%add the other textures 
    saveOptions.cropBorders = true;
    savePlot(fig, saveOptions);  

end

function [] = plotVisualResult(Ibase, Ioverlay, figTitle, labels, coordinates, cmap, hideColorbar, fig,saveOptions)
%% Visualization of malignancy score %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (nargin < 4)
        labels =[];
    end
    if (nargin < 5)
        coordinates = [];
    end
    if (nargin < 6)
        cmap = 'jet';
    end
    if (nargin < 7)
        hideColorbar = false;
    end 
    if (nargin < 8)
        fig = figure;
    else     
        figure(fig);
        clf(fig);
    end
    if (nargin < 9)
        saveOptions.SaveImage = false;
    end
    
    warning('off')

    alpha = 0.7;
    cmapSize = 100; % default size of 60 shows visible discretization
    if ischar(cmap)
        try
            cmap = eval([cmap '(' num2str(cmapSize) ');']);
        catch
            fprintf('Colormap ''%s'' is not supported. Using ''jet''.\n',cmapName);
            cmap = jet(cmapSize);
        end
    end
    colormap(cmap);

    clf(gcf);
    axes(gcf);
    
    imshow(Ibase)
    hold on 
    h = imagesc(Ioverlay);
    set(gca, 'CLim', [0 1]);
    colormap(cmap);
    hold off
    alphaI = Ioverlay > 0;
    set(h, 'AlphaData', alphaI * 0.5);
    if ~hideColorbar
        c = colorbar('location','southoutside', 'Ticks', [0 0.5 1], 'TickLabels', {'low', 'medium', 'high'});
        c.Label.String = 'Malignancy Probability';
        c.Label.FontSize = 10;
        c.Label.FontWeight = 'bold';
        c.Limits = [0,1];
        c.LimitsMode = 'manual';
        set(gcf,'Visible','on');
    end
        
    title(figTitle)
    if ~isempty(coordinates)
        hold on
        for i = 1:size(coordinates, 1)
            x = coordinates(i,1);
            y = coordinates(i,2);
            if (labels(i) == 1) 
                plot(x, y,'rx', 'MarkerSize', 10, 'MarkerEdgeColor', 'r', 'LineWidth', 3);
            else 
                plot(x, y,'go', 'MarkerSize', 10, 'MarkerEdgeColor', 'g', 'LineWidth', 3);
            end
        end
        hold off
    end
    
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 1, 1]);
    
    saveOptions.cropBorders = true;
    savePlot(fig, saveOptions);  
    warning('on')

end


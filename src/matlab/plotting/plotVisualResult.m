function [] = plotVisualResult(Ibase, Ioverlay, figTitle, labels, coordinates, cmap, hideColorbar, fig,saveOptions)
%% Visualization of malignancy score %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (nargin < 4)
        labels = {};
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
    
    key = {'Benign', 'Atypical', 'Malignant'};
    value = {'go', 'ms', 'r^'}; 
    lineStyleMap = containers.Map(key, value);
    
    value = {'g', 'm', 'r'}; 
    lineColorMap = containers.Map(key, value);

    warning('off')

    cmapSize = 100; % default size of 60 shows visible discretization
    if ischar(cmap)
        try
            cmap = eval([cmap '(' num2str(cmapSize) ');']);
        catch
            fprintf('Colormap ''%s'' is not supported. Using ''jet''.\n',cmapName);
            cmap = jet(cmapSize);
        end
    end

    clf(gcf);
    axes(gcf);
    
    imshow(Ibase)
    hold on 
    h = imagesc(Ioverlay);
    set(gca, 'CLim', [0 1]);
    colormap(cmap( 15:cmapSize, :));
    hold off
    alphaI = Ioverlay > 0;
    set(h, 'AlphaData', alphaI * 0.5);
    if ~hideColorbar
        c = colorbar('location','southoutside', 'Ticks', [0 0.5 1], 'TickLabels', {'low', 'medium', 'high'});
        c.Label.String = 'Malignancy Probability';
        c.Label.FontSize = 15;
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
            plot(x, y,lineStyleMap(labels{i}), 'MarkerSize', 12, 'MarkerEdgeColor', 'k', ...
                'MarkerFaceColor', lineColorMap(labels{i}),  'LineWidth', 3);
        end
        hold off
    end
    
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 1, 1]);
    
      
    hold on; 
    h = [];
    for i = 1:length(key)
        h(i) = plot(nan, nan, lineStyleMap(key{i}), 'MarkerSize', 12, 'MarkerEdgeColor', 'k', ...
                'MarkerFaceColor', lineColorMap(key{i}),  'LineWidth', 3, 'DisplayName', key{i});
    end
    hold off;
    
    hleg = legend(h);
    legPos = get(hleg,'position') ;
    set(hleg,'position', [legPos(1),legPos(2) - legPos(4), legPos(3),legPos(4) * 2.0]) ;
                
    saveOptions.cropBorders = true;
    savePlot(fig, saveOptions);  
    warning('on')

end


function [] = plotSegmentation(Iorig, segmentMask, Isegment, coordinates, fig,saveOptions)
%% Plot Segmentation results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (nargin < 5)
        fig = figure;
    else     
        figure(fig);
        clf(fig);
    end
    if (nargin < 6)
        saveOptions.SaveImage = false;
    end
    warning('off')

    h1 = subplot(1,2,1);
    clf(h1)
    imoverlay_medical(rgb2gray(Iorig), segmentMask, [], [], 'hot', [], h1);
    c = colorbar('location', 'westoutside','Ticks', linspace(0,1,5), 'TickLabels', linspace(0,100,5));
    c.Label.String = 'Channel Agreement (%)';
    c.Label.FontSize = 15;
    c.Label.FontWeight = 'bold';
    c.LimitsMode = 'manual';
    c.Limits = [0,1];    
    title('MSI Segments', 'FontSize', 15)

    h2 = subplot(1,2,2);
    colors = jet(size(coordinates,1));
    imshow(Isegment);
    hold on
    for i = 1:size(coordinates, 1)
        x = coordinates(i,1);
        y = coordinates(i,2);
        plot(x, y, '*', 'LineWidth', 2, 'MarkerSize', 7, 'MarkerFaceColor', colors(i,:));
    end
    hold off
    title('Final Segment', 'FontSize', 15);
    %plotName = strcat(plotName, '_segments');
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 1, 1]);

    saveOptions.cropBorders = true;
    savePlot(fig, saveOptions);
    warning('on')
        
end


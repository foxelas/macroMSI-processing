function [] = plotMSIWithPOI(Iorig, coordinates, fig)

colors = jet(size(coordinates, 1));
imshow(Iorig);
hold on
for i = 1:size(coordinates, 1)
    x = coordinates(i, 1);
    y = coordinates(i, 2);
    plot(x, y, '*', 'LineWidth', 2, 'MarkerSize', 7, 'MarkerFaceColor', colors(i, :));
end
hold off
title('Final Segment', 'FontSize', 15);
%plotName = strcat(plotName, '_segments');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 1, 1]);

savePlot(fig);
warning('on')

end
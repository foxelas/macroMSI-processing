function [] = plotPointsOnImage(baseImage, xPoints, yPoints, isCombination, fig)
%%PLOTPOINTSONIMAGE plots a set of point coordinates with text and markers
%%on a base image 
%
%   Usage: 
%   plotPointsOnImage(baseImage, xPoints, yPoints, true, fig);

if isCombination
    [yy, xx] = meshgrid(xPoints, yPoints);
    xx = xx(:);
    yy = yy(:);
else
    xx = xPoints; 
    yy = yPoints;
end
    
imshow(baseImage);
hold on;
for i = 1:length(xx)
    plot(xx(i), yy(i), 'rx', 'MarkerSize', 20, 'LineWidth', 5);
%%%%%%%%%%%%%%%%%%%%%%%%%%
%     textStr = sprintf('P%d(%d,%d)', i, xx(i), yy(i));
    textStr = sprintf('P%d', i);
    text(xx(i)- 50, yy(i) + 50, textStr);
end
hold off;

%set(gcf,  'units','normalized','outerposition',[0 0 1 1]);
if ~contains(getSetting('plotName'), 'point')
    setSetting('plotName', mkNewDir(getSetting('savedir'), getSetting('saveFolder'), 'pointsOnImage'));
end
savePlot(fig);

end 
function [] = plotYFromHSI(baseImage, figTitle, xPoints, yPoints, fig)

if ndims(baseImage) > 2
   imagesc(squeeze(baseImage(:,:,2)));
else
   imagesc(baseImage);
end 
if ~isempty(figTitle)
    title(figTitle)
end

if ~isempty(xPoints) && ~isempty(yPoints)
    [yy, xx] = meshgrid(xPoints,yPoints); 
    xx = xx(:);
    yy = yy(:);

    hold on;
    for i = 1:length(xx)
        plot(xx(i),yy(i),'rx', 'MarkerSize', 10, 'LineWidth', 5);
        textStr = sprintf('P%d at (%d,%d)', i,xx(i), yy(i));
        text(xx(i),yy(i), textStr);
    end 
    hold off;
    set(gcf,  'units','normalized','outerposition',[0 0 1 1]);
    plotName = getSetting('plotName');
    setSetting('plotName', strcat(plotName, '-points'));
end 

colorbar; 
savePlot(fig);

end


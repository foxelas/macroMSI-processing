function [] = plotGFCHistogram(gfcs, fig,saveOptions)

    
if (nargin < 2)
    fig = figure;
else     
    figure(fig);
    clf(fig);
end
if (nargin < 3)
    saveOptions.SaveImage = false;
end
    
histogram(gfcs, 7);
ax  = get(gca);
ax.FontSize = 18;
xlabel('Goodness-Of-Fit Criterion', 'FontSize', 20);
ylabel('Number of POIs', 'FontSize', 20);
ylim([0, 80]);

savePlot(fig, saveOptions);

end
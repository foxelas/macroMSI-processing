function [] = plotGFCHistogram(gfcs, fig)
%     PLOTGFCHISTOGRAM plots the GFC histogram
%
%     Usage:
%     plotGFCHistogram(gfcs, fig)

histogram(gfcs, 7);
ax = get(gca);
ax.FontSize = 18;
xlabel('Goodness-Of-Fit Criterion', 'FontSize', 20);
ylabel('Number of POIs', 'FontSize', 20);
ylim([0, 80]);

savePlot(fig);

end
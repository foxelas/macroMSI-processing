function [] = plotNMSE2(x1, x2, y, x1lab, x2lab, ylab, figTitle, fig)

warning('off');

y = reshape(y, [length(x1), length(x2)]);
[X1, X2] = meshgrid(x1, x2);
s = surf(X1, X2, y', 'FaceAlpha', 0.5);
s.Parent.View = [-63.0400, 25.3831];
s.EdgeColor = 'interp';
ax = gca;
ax.FontSize = 20;

xlabel(x1lab, 'FontSize', 30, 'Interpreter', 'latex');
xticks(x1);
ylabel(x2lab, 'FontSize', 30, 'Interpreter', 'latex');
yticks(x2);
zlabel(ylab, 'FontSize', 30);
%title(figTitle, 'FontSize', 15);
c = colorbar('Location', 'northoutside');
c.LimitsMode = 'manual';
c.Limits = [0, 0.1];
pause;
savePlot(fig);

warning('on');
end
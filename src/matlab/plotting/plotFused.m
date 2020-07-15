function [] = plotFused(left, right, pointsLeft, pointsRight, figTitle, fig)

warning('off');

fig = figure(fig);
ax = axes;
showMatchedFeatures(left,right,pointsLeft,pointsRight, 'Parent', ax);
title(figTitle);
legend('ptsFixed','ptsMoving');
pause(0.1)
savePlot(fig);

warning('on');
end

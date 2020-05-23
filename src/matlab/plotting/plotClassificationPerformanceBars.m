function [] = plotClassificationPerformanceBars(performance, barLabels, legends, labely, lims, fig)

%% Plot bar performance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin < 5)
    lims = [70, 100];
end

color = hsv(3);
%     if (saveInBW)
%         color(1,:)=HInt2RGB(1,100); % red, darkest
%         color(2,:)=HInt2RGB(3,64); % green, less dark
%         color(3,:)=HInt2RGB(7,10); % blue cyan, lightest
%     else
%         color = hsv(3);
%         %color = flip(color);
%     end
if any(contains(legends, 'SVM')) && any(contains(barLabels, 'LBP'))
    legTitle = 'Classifier';
    %figTitleCur = 'Cross Validated AUC';
    labelx = 'Feature Vector';

elseif any(contains(legends, 'fixed')) && any(contains(barLabels, 'LBP'))
    legTitle = 'Tissue';
    %figTitleCur = 'Cross Validated AUC';
    labelx = 'Feature Vector';
else
    legTitle = 'Tissue';
    labelx = 'Classifier';
    labely = strcat([figTitleCur, '(%)']);
    %figTitleCur = strcat(['Cross Validated ', figTitleCur]);
end
hasLegend = true;
plotFunWrapper( fig, @plotBars, barLabels, performance*100, color, '',...
    labelx, labely, legends, legTitle, lims, 15, hasLegend);
savePlot(fig);

end

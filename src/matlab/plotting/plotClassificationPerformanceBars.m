function [] = plotClassificationPerformanceBars(performance,lineNames,figTitle,fig,saveOptions)
%% Plot bar performance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if (nargin < 4)
        fig = figure;
    else 
        figure(fig);
        clf(fig);
    end

    if (nargin < 5)
        saveOptions.SaveImage = false;
    end
    
    if (saveInBW)
        color(1,:)=HInt2RGB(1,100); % red, darkest 
        color(2,:)=HInt2RGB(3,64); % green, less dark 
        color(3,:)=HInt2RGB(7,10); % blue cyan, lightest
    else 
        color = hsv(3);
        %color = flip(color);
    end

    plotN = length(figTitle)-1;
    superTitle = figTitle{1};
    figTitle = figTitle(2:end); 
    for ploti=1:plotN 
        subplot(plotN,1,ploti);
        figTitleCur = figTitle{ploti};
        performanceCur = performance{ploti};
        hasLegend = ploti == 1;
        if contains(figTitleCur, 'classifier')
            legTitle = 'Classifier';
            figTitleCur = 'Cross Validated AUC';
            labelx = 'Feature Vector';
            labely = 'ROC AUC';
            cs = {'SVM', 'KNN', 'RF'};
            lims = [70, 100];
        elseif contains(figTitleCur, 'fixing')
            legTitle = 'Tissue';
            figTitleCur = 'Cross Validated AUC';
            labelx = 'Feature Vector';
            labely = 'ROC AUC';
            cs = {'Unfixed', 'Fixed', 'Mixed'};
            lims = [70, 100];
        else 
            legTitle = 'Tissue';
            cs = {'Unfixed', 'Fixed', 'Mixed'};
            labelx = 'Classifier';
            labely = strcat([figTitleCur, '(%)']);
            figTitleCur = '';
            %figTitleCur = strcat(['Cross Validated ', figTitleCur]);
            lims = [50, 100];
        end
        barPlot = GetBarPlot(lineNames, performanceCur * 100, color, figTitleCur, labelx, labely, cs, legTitle, lims, 9, hasLegend);
    end 
    suptitle(superTitle);
    saveOptions.plotName = fullfile(saveOptions.savedir, superTitle);
      
    savePlot(fig, saveOptions);

end

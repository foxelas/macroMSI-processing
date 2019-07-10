function [] = plotROCs(foldPerformance, averagePerformance, figTitle, fig,saveOptions)

    if (nargin < 4)
        fig = figure;
    else     
        figure(fig);
        clf(fig);
    end
    if (nargin < 5)
        saveOptions.SaveImage = false;
    end
    
    hold on
    folds = length(foldPerformance);
    c = colormap(folds + 1);% colormap(lines);    
    for i = 1:folds*(folds < 11)
        if ~isempty(foldPerformance(i).AUC)
            plot(foldPerformance(i).ROCX, foldPerformance(i).ROCY, 'DisplayName', sprintf('Fold %d (AUC = %.3f)', i, foldPerformance(i).AUC), 'LineWidth', 2);
        end
    end
    plot(0:0.1:1, 0:0.1:1, 'k--', 'DisplayName', 'Chance');
    if ~isempty(averagePerformance.ROCX)
        plot(averagePerformance.ROCX(:,1), averagePerformance.ROCY(:,1), 'm-*',  'DisplayName', sprintf('Average (AUC = %.3f)', averagePerformance.AUC(1)), 'LineWidth', 2);
%         shadedErrorBar(performance.ROCX(:,1), performance.ROCY(:,1), abs(performance.ROCY(:,1) - performance.ROCY(:,2:3)),'lineprops','m-*');
    end
    hold off 
    xlim([-0.1, 1.1]);
    ylim([-0.1, 1.1]);
    xlabel('False positive rate');
    ylabel('True positive rate');
    h = findobj(gca,'Type','line');
    hh = h(contains({h.DisplayName}, {'AUC', 'Chance'}));
    legend(hh, 'Location', 'eastoutside');
    title(figTitle);
    set(gcf, 'Position', get(0, 'Screensize'));
    
    savePlot(fig, saveOptions);
end


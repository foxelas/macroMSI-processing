function [barPlot] = GetBarPlot(lineNames, values, color, figTitle, labelx, labely, categories, legendTitle, ylims, ticklabelsize, hasLegend)

            barPlot = bar(categorical(lineNames), values, 'FaceColor', 'flat');
            for k = 1:length(categories)
                barPlot(k).CData = color(k,:);
            end
            hleg = legend(categories, 'FontSize', 9.5, 'Location', 'eastoutside');
            htitle = get(hleg,'Title');
            set(htitle,'String',legendTitle)
            if ~(hasLegend)
                hleg.Visible = 'off';
            end
            ax = gca;
            ax.FontSize = ticklabelsize; 
            ylim(ylims);
            xlabel(labelx, 'FontSize', 12);
            ylabel(labely, 'FontSize', 12);
            if ~strcmp(figTitle, '')
                title(figTitle, 'FontSize', 16); % should be after gca, set to work
            end
end
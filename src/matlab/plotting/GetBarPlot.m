function [barPlot] = GetBarPlot(lineNames, values, color, figTitle, labelx, labely, categories, legendTitle, ylims, ticklabelsize, hasLegend, isStacked)

            if (nargin < 12)
                isStacked = false;
            end
            
            if (length(lineNames) ~= size(values, 1))
                values  = values';
            end
          
            if (isStacked) 
                barPlot = bar(categorical(lineNames), values, 'stacked', 'FaceColor', 'flat');
            else 
                barPlot = bar(categorical(lineNames), values, 'FaceColor', 'flat');
            end
            for k = 1:length(categories)
                barPlot(k).CData = color(k,:);
            end
            hleg = legend(categories, 'FontSize', 20, 'Location', 'eastoutside');
            htitle = get(hleg,'Title');
            set(htitle,'String',legendTitle)
            if ~(hasLegend)
                hleg.Visible = 'off';
            end
            ax = gca;
            ax.FontSize = ticklabelsize; 
            ylim(ylims);
            xlabel(labelx, 'FontSize', 20);
            ylabel(labely, 'FontSize', 20);
            if ~strcmp(figTitle, '')
                title(figTitle, 'FontSize', 16); % should be after gca, set to work
            end
end
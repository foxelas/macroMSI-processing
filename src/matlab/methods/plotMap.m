function plotMap(I, mask, cmap, hideColorbar, barTitle, limits, ticks, tickLabels, figTitle, fig)
% PLOTMAP plots a map with a colorbar (optional)
%
%     Usage: 
%     plotMap(msi, 'jet', false, 'Scale', [0,1], [0, 0.5, 1], {'low', 'medium', 'high'}, 'Map', 1);

if isempty(cmap)
    cmap = 'jet';
end 

if isempty(hideColorbar)
    hideColorbar = false; 
end 

cmapSize = 100; % default size of 60 shows visible discretization
if ischar(cmap)
    try
        cmap = eval([cmap, '(', num2str(cmapSize), ');']);
    catch
        fprintf('Colormap ''%s'' is not supported. Using ''jet''.\n', cmapName);
        cmap = jet(cmapSize);
    end
end

clf(gcf);
axes(gcf);

warning('off')

h = imagesc(I);
axis off;
alphaI = mask > 0;
opacity = 1;
set(h, 'AlphaData', alphaI*opacity);
if isempty(limits)
    importantPixels = I(mask & abs(I) < Inf);
    limits = [min(importantPixels(:)), max(importantPixels(:))];
end
if ~isempty(limits) && ~hideColorbar
    set(gca, 'CLim', limits);
end

colormap(cmap(15:cmapSize, :)); 

if ~hideColorbar
    c = colorbar('location', 'southoutside');
    c.Label.String = barTitle;
    c.Label.FontSize = 15;
    c.Label.FontWeight = 'bold';
    if ~isempty(limits)
        c.Limits = limits ; % [0, 1]
        c.LimitsMode = 'manual';
    end
    if ~isempty(ticks) 
        c.Ticks = ticks; % [0, 0.5, 1]
    end 
    if ~isempty(tickLabels)
        c.TickLabels = tickLabels; %  {'low', 'medium', 'high'}
    end 
    set(gcf, 'Visible', 'on');
end

if ~isempty(figTitle)
    title(figTitle, 'FontSize', 15);
end

savePlot(fig);
warning('on')

end 
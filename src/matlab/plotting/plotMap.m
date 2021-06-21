function plotMap(I, mask, cmap, hideColorbar, barTitle, limits, figTitle, fig)
% PLOTMAP plots a map with a colorbar (optional)
%
%     Usage:
%     plotMap(msi, 'jet', false, 'Scale', [0,1], 'Map', 1);

if (nargin < 7)
    figTitle = '';
    fig = [];
end

hasNumberedScale = getSetting('hasNumberedScale');
if hasNumberedScale 
    tickLabels =  {'low', 'medium', 'high'};
    ticks = [0, 0.5, 1];
else
    tickLabels = [];
    ticks = [];
end 

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

% clf(gcf);
% axes(gcf);
axes(gca);

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
    if ~isempty(barTitle)
        c.Label.String = barTitle;
        c.Label.FontSize = 15;
        c.Label.FontWeight = 'bold';
    end

    if ~isempty(limits)
        c.Limits = limits; % [0, 1]
        c.LimitsMode = 'manual';
    end
    if ~isempty(ticks)
        c.Ticks = ticks; % [0, 0.5, 1]
    end
    if ~isempty(tickLabels)
        c.TickLabels = tickLabels; %  {'low', 'medium', 'high'}
    end
    c.Ticks = []; 
    
    set(gcf, 'Visible', 'on');
end

if ~isempty(figTitle)
    title(figTitle, 'FontSize', 15);
end

savePlot(fig);
warning('on')

end
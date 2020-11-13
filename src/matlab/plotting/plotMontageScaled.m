function[] = plotMontageScaled(imageList, imageMasks, imageNames, filename, fig)

% Doesn't work well with many rows (e.g. above 4)

n = numel(imageList);
rows = floor(n/2);
cols = 2;
mask1 = imageMasks{1};
mask2 = imageMasks{2};

warning('off')
saveImages = getSetting('saveImages');
setSetting('saveImages', false);
cropBorders = getSetting('cropBorders');
setSetting('cropBorders', false);

set(gcf, 'Position', [500, 50, 760, 700]);
t = tiledlayout(rows, cols);
for i = 1:rows
    curIdx = i * cols - 1;
    map1 = imageList{curIdx};
    map2 = imageList{curIdx+1};
    mapLimits = []; % [min(min(map1(:)), min(map2(:))), max(max(map1(:)), max(map2(:)))];

    h(1) = nexttile; 
    plotMap(map1, mask1, [], true, [], mapLimits); %imageNames{curIdx}
    if i == 1
        title('Unfixed');
    end 
    if i == rows 
        pos = get(h(1), 'OuterPosition');
        text(pos(1) - pos(3)/2, -0.22 , 'Relative Chromophore Distribution', 'FontSize',15, 'Units', 'normalized'); 
    end 
    h(2) = nexttile; 
    plotMap(map2, mask2, [], true, [], mapLimits); % imageNames{curIdx + 1}
    if i == 1
        title('Fixed');
    end
    if i == rows 
        cb = colorbar(h(2), 'Location', 'southoutside'); 
    end 
end

t.Padding = 'compact';
t.TileSpacing = 'none';

setSetting('saveImages', saveImages);
savedir = getSetting('savedir');
mapdir = getSetting('map');
setSetting('plotName', fullfile(savedir, mapdir, strcat(filename, '_', 'montage.png')));
savePlot(fig);
setSetting('cropBorders', cropBorders);

warning('on');

end
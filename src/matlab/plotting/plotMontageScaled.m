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

set(gcf, 'Position', [407, -54, 675, 738]);
t = tiledlayout(rows, cols);
for i = 1:rows
    curIdx = i * cols - 1;
    map1 = imageList{curIdx};
    map2 = imageList{curIdx+1};
    mapLimits = [min(min(map1(:)), min(map2(:))), max(max(map1(:)), max(map2(:)))];

    nexttile; %subplot(rows, cols, curIdx);
    plotMap(map1, mask1, [], false, [], mapLimits); %imageNames{curIdx}
    nexttile; %subplot(rows, cols, curIdx + 1);
    plotMap(map2, mask2, [], false, [], mapLimits); % imageNames{curIdx + 1}
end
t.Padding = 'none';
t.TileSpacing = 'none';

setSetting('saveImages', saveImages);
savedir = getSetting('savedir');
mapdir = getSetting('map');
setSetting('plotName', fullfile(savedir, mapdir, strcat(filename, '_', 'montage.png')));
savePlot(fig);
setSetting('cropBorders', cropBorders);

warning('on');

end
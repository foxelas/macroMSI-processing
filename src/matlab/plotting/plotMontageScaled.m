function[] = plotMontageScaled(imageList, imageMasks, filename, fig)

% Doesn't work well with many rows (e.g. above 4)

n = numel(imageList);
rows = floor(n/2);
cols = 2;
mask1 = imageMasks{1};
mask2 = imageMasks{2};

nameparts = strsplit(filename, '_');
sampleName = nameparts{1};
if strcmp(nameparts{3}, 'Mel')
    text1 = 'Relative Melanin Concentration';
elseif strcmp(nameparts{3}, 'Hb')
    text1 = 'Relative Hemoglobin Concentration';
else 
    text1 = '';
end

if strcmp(nameparts{2}, 'Hb')
    text2 = strcat(sampleName, ' HHbT ROI');
elseif strcmp(nameparts{2}, 'Mel')
    text2 = strcat(sampleName, ' HM ROI');
elseif strcmp(nameparts{2}, 'Norm')
    text2 = strcat(sampleName, ' Norm ROI');
else 
    text2 = '';
end 

warning('off')
saveImages = getSetting('saveImages');
setSetting('saveImages', false);
cropBorders = getSetting('cropBorders');
setSetting('cropBorders', false);

hasLimits = getSetting('hasLimits');
% hasNumberedScale = getSetting('hasNumberedScale');

set(gcf, 'Position', [500, 50, 760, 700]);
t = tiledlayout(rows, cols);
for i = 1:rows
    curIdx = i * cols - 1;
    map1 = imageList{curIdx};
    map2 = imageList{curIdx+1};
    if hasLimits
        mapLimits = [0,1];
%         mapLimits = [min(min(map1(:)), min(map2(:))), max(max(map1(:)), max(map2(:)))];
    else 
        mapLimits = []; 
    end
    
    h(1) = nexttile; 
    plotMap(map1, mask1, [], true, [], mapLimits); %imageNames{curIdx}
    if i == 1
        title('Unfixed');
    end 
    if i == rows 
        pos = get(h(1), 'OuterPosition');
        text(pos(1) - pos(3)/2, -0.22 , text1, 'FontSize',15, 'Units', 'normalized'); 
        text(pos(1) - pos(3)/2, -0.32 , text2, 'FontSize',15, 'Units', 'normalized'); 
    end 
    h(2) = nexttile; 
    plotMap(map2, mask2, [], true, [], mapLimits); % imageNames{curIdx + 1}
    if i == 1
        title('Fixed');
    end
    if i == rows 
%         figHandles = findobj('Type', 'figure');
%         figHandle = figHandles(find([figHandles.Number] == 2));
        caxis([0, 1]);
        cb = colorbar( 'Location', 'southoutside'); 
        cb.LimitsMode = 'manual';
        cb.Limits =  [0, 1]; %mapLimits;
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
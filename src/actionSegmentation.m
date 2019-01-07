defaultfc = [450, 465, 505, 525, 575, 605, 630];
% Colors defined from https://academo.org/demos/wavelength-to-colour-relationship/
bandColors = [  0,70,255   ; 
                0,146,255  ;
                0,255,84   ;
                74,255,0   ;
                240,255,0  ;
                255,173,0  ;
                255,79,0   ;
                255,255,255
              ];
bandColors = bandColors./255;

[~, idx] = unique({ID.IMG});

for k = idx'
    
    % Give a name
    [options.saveOptions.plotName, sampleName] = generateName(options, 'plot+save', data(ID(k).Representative), ID(k));
    
    % Retrieve whole MSI
    MSI = readMSI({data(ID(k).Data).File}, []); %     I = squeeze(MSI(1,:,:,:));
    
    g = permute(valueSelect(MSI, 'adjusted'), [2, 3, 1]);
    
    x = ID(k).Originx;
    y = ID(k).Originy;
    
    fig1 = figure(1);
    imshow(squeeze(MSI(4,:,:,:)));
    hold on
    plot(x, y, 'r*', 'LineWidth', 2, 'MarkerSize', 15);
    hold off
    
    mm = regionGrowing2(g, [y, x, 4]);


    [m, n, ~] = size(g);
    mask = zeros(m, n, 3);
    for i = 1:7  
        I2 = squeeze(g(:,:,i));
        [~, maskTmp]= regionGrowing(I2, [y, x]);
        mask = mask + cat(3, bandColors(i,1) * maskTmp, bandColors(i,2) * maskTmp, bandColors(i,3) * maskTmp);
    end
        
    fig2 = figure(2);
    colormap(bandColors);
    rgbImage = cat(3, I2, I2, I2);
    masked = I2 + mask;
    imshow(masked);
    c = colorbar('location','southoutside', 'Ticks', linspace(0,1,9),...
         'TickLabels',{'450','465','505','525','575','605', '630', 'All', ''});
    c.Label.String = 'Respective MSI band (nm)';
    
    if (options.saveOptions.saveImages)
        set(0, 'CurrentFigure', fig1)
        export_fig([options.saveOptions.plotName, '_origin.jpg'] , '-jpg')
        set(0, 'CurrentFigure', fig2)
        export_fig([options.saveOptions.plotName, '.jpg'] , '-jpg','-native')
    end
end
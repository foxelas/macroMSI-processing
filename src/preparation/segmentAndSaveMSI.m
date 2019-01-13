function [] = segmentAndSaveMSI(options, ID, data)

tic;

whiteMSIStruct = struct('Name', {}, 'Index', [], 'MSI', []);
darkMSIStruct = struct('Name', {}, 'Index', [], 'MSI', []);
MSIStruct = struct('Name', {}, 'Index', [], 'MSI', [], 'Mask', []);
appendingIdx = 1;

for k = 1:length(ID)
    % Give a name
    [options.saveOptions.plotName, ~] = generateName(options, 'plot+save', data(ID(k).Representative), ID(k));
    name = generateName([], 'image', data(ID(k).Representative), ID(k));

    files = {data(ID(k).Data).File};
    x = ID(k).Originx;
    y = ID(k).Originy;
    
    [MSI, patchMask, whiteReference, darkReference] = SegmentAndSave(files, x, y, options);
    
    if ~(isempty(whiteReference))
        whiteMSIStruct(appendingIdx) = struct('Name', name, 'Index', appendingIdx, 'MSI', whiteReference);
    end
    
    if ~(isempty(darkReference))
        darkMSIStruct(appendingIdx) = struct('Name', name, 'Index', appendingIdx, 'MSI', darkReference);
    end
    
    if ~(isempty(MSI))
        MSIStruct(appendingIdx) = struct('Name', name, 'Index', appendingIdx, 'MSI', MSI, 'Mask', patchMask);
        appendingIdx = length(MSIStruct) + 1;
    end
      
end
t = toc;
fprintf('Action elapsed %.2f hours.\n', t / (60 * 60));

m = matfile(generateName(options, 'matfilein'), 'Writable', true);
m.MSIStruct = orderfields(MSIStruct);
m.WhiteMSIStruct = orderfields(whiteMSIStruct);
m.DarkMSIStruct = orderfields(darkMSIStruct);

fprintf('Finished reading all MSI data according to ID file (Region Growing case). Saving in %s.\n', options.systemdir);
end 

function [patch, patchMask, whiteI, darkI] = SegmentAndSave( files, x, y, options, accTheta, regionRadius)
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
    bands = length(defaultfc);
    bandColors = bandColors./255;
    bandWeight = 1 / bands;
    
    if (nargin < 5) || isempty(accTheta)
        accTheta = 0.6;
    end
    if (nargin < 6) || isempty(regionRadius)
        regionRadius = 50;
    end

    % Retrieve whole MSI
    [MSI, whiteI, darkI] = readMSI(files, []); %     I = squeeze(MSI(1,:,:,:));
    
    g = permute(valueSelect(MSI, 'adjusted'), [2, 3, 1]);
    
    [m, n, ~] = size(g);
    maskForFig = zeros(m, n, 3); % RGB colored mask to show which MSI band resulted in which region
    mask = zeros(m, n); % Final mask for the region 
  
    for i = 1:bands  
        I2 = squeeze(g(:,:,i));
        [~, maskTmp]= regionGrowing(I2, [y, x], [], regionRadius);
        maskForFig = maskForFig + cat(3, bandColors(i,1) * maskTmp, bandColors(i,2) * maskTmp, bandColors(i,3) * maskTmp);
        mask = mask + bandWeight * maskTmp; 
    end
    
    mask = mask > accTheta; 
    if sum(mask(:)) < 9 %if the region is too small, add the neigborhood of the region seed pixel
        patchX = (x-2):(x+2);
        patchY = (y-2):(y+2);
        mask(patchY, patchX) = 1;
    end
    [r, c] = find(mask);
    patchY = min(r):max(r);
    patchX = min(c):max(c);
    patch = MSI(:, patchY, patchX, :);
    patchMask = mask(patchY, patchX);
    
    if false %(options.showImages)
        fig1 = figure(1);
        imshow(whiteI + mask); %imshow(squeeze(MSI(4,:,:,:)));
        hold on
        plot(x, y, 'r*', 'LineWidth', 2, 'MarkerSize', 15);
        hold off
    
        fig2 = figure(2);
        colormap(bandColors);
%         rgbImage = cat(3, I2, I2, I2);
        masked = I2 + maskForFig;
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
    
    if ~isempty(whiteI)
        whiteI = whiteI(patchY, patchX,:);
    end
    if ~isempty(darkI)
        darkI = darkI(patchY, patchX,:);
    end
end

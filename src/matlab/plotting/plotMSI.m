function [] = plotMSI(raw,fig,saveOptions)

    if (nargin < 2)
        fig = figure;
    else 
        figure(fig);
        clf(fig);
    end

    if (nargin < 3)
        saveOptions.SaveImage = false;
    end 

    warning('off');
    saveOptions.plotName  = fullfile(saveOptions.savedir, 'MSI', 'msi');

    msi = raw2msi(raw, 'extended');    
    msi = permute(msi, [2, 3, 1]);
    montage(msi, 'Size', [3, 3]);
    saveOptions.cropBorders = true;
    pause(0.1)
    savePlot(fig, saveOptions);
    
    origPlotName = saveOptions.plotName;
    for i=1:size(msi, 3)
        imshow(squeeze(msi(:,:,i)));
        saveOptions.plotName = strcat(origPlotName,'_', num2str(i));
        saveOptions.cropBorders = true;
        pause(0.1)
        savePlot(fig, saveOptions);
    end
    

    
    warning('on');
end

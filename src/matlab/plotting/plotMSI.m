function [] = plotMSI(raw,fig,saveOptions, isSeparateImages)

    if (nargin < 2)
        fig = figure;
    else 
        figure(fig);
        clf(fig);
    end

    if (nargin < 3)
        saveOptions.SaveImage = false;
    end 

    if (nargin < 4)
        isSeparateImages = false;
    end 
    
    warning('off');
    saveOptions.plotName  = fullfile(saveOptions.savedir, 'MSI', 'msi');
    saveOptions.cropBorders = true;

    if ndims(raw) > 3 
        msi = raw2msi(raw, 'extended');    
        %msi = permute(msi, [2, 3, 1]);
    else 
        msi = raw; 
    end 
    channels = size(msi, 1);

    if (isSeparateImages)
        origPlotName = saveOptions.plotName;
        for i=1:channels
            imshow(squeeze(msi(i,:,:)));
            saveOptions.plotName = strcat(origPlotName,'_', num2str(i));
            pause(0.1)
            savePlot(fig, saveOptions);
        end
    else
        imageList = num2cell(msi, [2, 3]);
        imageList = cellfun(@squeeze, imageList, 'un', 0);
        montage(imageList, 'Size', [2, ceil(channels/2)]);
        pause(0.1)
        savePlot(fig, saveOptions);
    end 

    warning('on');
end

function [] = plotMSI(raw, fig, isSeparateImages)

if (nargin < 3)
    isSeparateImages = false;
end

warning('off');
setSetting('cropBorders', true);
plotName = fullfile(getSetting('savedir'), 'MSI', 'msi');
setSetting('plotName', plotName);

if ndims(raw) > 3
    msi = raw2msi(raw, 'extended');
    %msi = permute(msi, [2, 3, 1]);
else
    msi = raw;
end
channels = size(msi, 1);

if (isSeparateImages)
    origPlotName = plotName;
    for i = 1:channels
        imshow(squeeze(msi(i, :, :)));
        setSetting('plotName', strcat(origPlotName, '_', num2str(i)));
        pause(0.1)
        savePlot(fig);
    end
else
    imageList = num2cell(msi, [2, 3]);
    imageList = cellfun(@squeeze, imageList, 'un', 0);
    montage(imageList, 'Size', [2, ceil(channels/2)]);
    pause(0.1)
    savePlot(fig);
end

warning('on');
end

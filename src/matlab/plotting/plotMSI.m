function [] = plotMSI(raw, isSeparateImages, fig)

if isempty(isSeparateImages)
    isSeparateImages = false;
end

warning('off');
setSetting('cropBorders', true);
if ndims(raw) > 3
    msi = raw2msi(raw, 'adjusted');
    %msi = permute(msi, [2, 3, 1]);
else
    msi = raw;
end
channels = size(msi, 1);

if (isSeparateImages)
    plotName = fullfile(getSetting('savedir'), 'MSI', 'msi');
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
    if ~exist('fig', 'var')
        fig = gcf; 
    end 
    savePlot(fig);
end

warning('on');
end

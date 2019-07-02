function [] = visualizeLBP(raw,whiteReference,specimenMask, g, saveOptions)

    msi = raw2msi(raw, 'extended');
    [~,M,N] = size(msi);

    sumLBPImage = zeros(M-2, N-2);
    catLBPImage = [];
    for i = 1:msibands
        msi(i,:,:) = squeeze(msi(i,:,:)) .* specimenMask;
        lbps = lbp(squeeze(msi(i,:,:)),scale,neighbors,mapping, 'e');
        sumLBPImage = sumLBPImage + im2double(lbps);
        lbps = im2double(lbps ./10);
        catLBPImage = [catLBPImage ; lbps];
    end

    mmLBPImage = im2double(lbp(msi, scale, neighbors, mapping, 'e') ./ 40);
    LBPImage = im2double(lbp(rgb2gray(whiteReference) .* specimenMask,scale,neighbors,mapping, 'e'));

    type = 'SumLBP';
    saveOptions.plotName = fullfile(saveOptions.savedir, '7-LBPVisualisation', strcat( 'lbp_', num2str(g), '_', type));
    plotLBP(sumLBPImage, type, 1, saveOptions);
    
    type = 'CatLBP';
    saveOptions.plotName = fullfile(saveOptions.savedir, '7-LBPVisualisation', strcat( 'lbp_', num2str(g), '_', type));
    plotLBP(catLBPImage .* 100, type, 2, saveOptions);
    
    type = 'MMLBP';
    saveOptions.plotName = fullfile(saveOptions.savedir, '7-LBPVisualisation', strcat( 'lbp_', num2str(g), '_', type));
    plotLBP(mmLBPImage .* 100, type, 3, saveOptions);
    
    type = 'LBP';
    saveOptions.plotName = fullfile(saveOptions.savedir, '7-LBPVisualisation', strcat( 'lbp_', num2str(g), '_', type));
    plotLBP(LBPImage, type, 4, saveOptions);

end 
function [] = visualizeLBP(raw, whiteReference, specimenMask, g, saveOptions)

msi = raw2msi(raw, 'extended');
[B, M, N] = size(msi);
neighbors = 8;
mapping = getmapping(neighbors, 'riu2');
outputFolderMap = getOutputDirectoryMap();

figure(1);
figure(2);
set(gcf, 'Position', get(0, 'Screensize'));
figure(3);
set(gcf, 'Position', get(0, 'Screensize'));
figure(4);

for scale = 1:3
    catLBPImage = [];
    sumLBPImage = zeros(M-2*scale, N-2*scale);
    for i = 1:B
        msi(i, :, :) = squeeze(msi(i, :, :)) .* specimenMask;
        lbps = lbp(squeeze(msi(i, :, :)), scale, neighbors, mapping, 'e');
        sumLBPImage = sumLBPImage + im2double(lbps);
        lbps = im2double(lbps./10);
        catLBPImage = [catLBPImage; lbps];
    end


    mmLBPImage = im2double(lbp(msi, scale, neighbors, mapping, 'e')./40);
    LBPImage = im2double(lbp(rgb2gray(whiteReference).*specimenMask, scale, neighbors, mapping, 'e'));

    type = 'SumLBP';
    saveOptions.plotName = fullfile(saveOptions.savedir, outputFolderMap('lbpVisualization'), strcat('lbp_', num2str(g), '_', type, '_', num2str(scale)));
    plotLBP(sumLBPImage, type, 1, saveOptions);

    type = 'CatLBP';
    saveOptions.plotName = fullfile(saveOptions.savedir, outputFolderMap('lbpVisualization'), strcat('lbp_', num2str(g), '_', type, '_', num2str(scale)));
    plotLBP(catLBPImage.*100, type, 2, saveOptions);

    type = 'MMLBP';
    saveOptions.plotName = fullfile(saveOptions.savedir, outputFolderMap('lbpVisualization'), strcat('lbp_', num2str(g), '_', type, '_', num2str(scale)));
    plotLBP(mmLBPImage.*100, type, 3, saveOptions);

    type = 'LBP';
    saveOptions.plotName = fullfile(saveOptions.savedir, outputFolderMap('lbpVisualization'), strcat('lbp_', num2str(g), '_', type, '_', num2str(scale)));
    plotLBP(LBPImage, type, 4, saveOptions);


end
end
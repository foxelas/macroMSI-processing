function [] = visualizeLBP(raw, whiteReference, specimenMask, g)
%     VISUALIZELBP visualizes the lbp features
%
%     Usage:
%     visualizeLBP(raw, whiteReference, specimenMask, g)

msi = raw2msi(raw, 'extended');
[B, M, N] = size(msi);
neighbors = 8;
mapping = getmapping(neighbors, 'riu2');

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
    setSetting('plotName', fullfile(getSetting('savedir'), getSetting('lbpVisualization'), strcat('lbp_', num2str(g), '_', type, '_', num2str(scale))));
    plots(1, @plotLBP, sumLBPImage, type);

    type = 'CatLBP';
    setSetting('plotName', fullfile(getSetting('savedir'), getSetting('lbpVisualization'), strcat('lbp_', num2str(g), '_', type, '_', num2str(scale))));
    plots(2, @plotLBP, catLBPImage.*100, type);

    type = 'MMLBP';
    setSetting('plotName', fullfile(getSetting('savedir'), getSetting('lbpVisualization'), strcat('lbp_', num2str(g), '_', type, '_', num2str(scale))));
    plots(3, @plotLBP, mmLBPImage.*100, type);

    type = 'LBP';
    setSetting('plotName', fullfile(getSetting('savedir'), getSetting('lbpVisualization'), strcat('lbp_', num2str(g), '_', type, '_', num2str(scale))));
    plots(4, @plotLBP, LBPImage, type);


end
end
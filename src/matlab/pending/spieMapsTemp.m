makeExcel = false;
if makeExcel 
    clear all ;

    [meanMaps, Ncc2d, kld2d, hi2d, emd2d ]= getAccumArrays( @to2dMatrixMean);
    [stdMaps, Ncc2dS, kld2dS, hi2dS, emd2dS ]= getAccumArrays( @to2dMatrixStd);

    [imgArray, nccArray, kldArray, hiArray, emdArray] = formatArrays(true, meanMaps, stdMaps, Ncc2d, kld2d, hi2d, emd2d, Ncc2dS, kld2dS, hi2dS, emd2dS);

    writecell(imgArray,'roiArrays.xlsx','Sheet',1,'Range','A1');
    writecell(nccArray,'roiArrays.xlsx','Sheet',1,'Range','A10');
    writecell(kldArray,'roiArrays.xlsx','Sheet',1,'Range','A20');
    writecell(hiArray,'roiArrays.xlsx','Sheet',1,'Range','A30');
    writecell(emdArray,'roiArrays.xlsx','Sheet',1,'Range','A40');

    save('roiArrays.mat');

    fileID = fopen('tables.txt','w');
    fprintf(fileID,table2latex(imgArray), '\n\n');
    fprintf(fileID,table2latex(nccArray), '\n\n');
    fprintf(fileID,table2latex(kldArray), '\n\n');
    fprintf(fileID,table2latex(hiArray), '\n\n');
    fprintf(fileID,table2latex(emdArray), '\n\n');
    fclose(fileID);

end 

makePlots = true; 
if makePlots 
    plotGraphs();
end 

function [yMeans] = getMeans(y)
    yMeans = cell2mat(cellfun(@(x) mean(x, 'all', 'omitnan'), y, 'UniformOutput', false));
end

function [mm] = to2dMatrixMean(yy)
    mm = mean(zeros2Nan(yy), 3, 'omitnan');
end 

function [mm] = to2dMatrixStd(yy)
    mm = std(zeros2Nan(yy), 0, 3, 'omitnan');
end 

function [nanyy] = zeros2Nan(yy)
    nanyy = yy; 
    for i =1:size(yy,1)
        for j=1:size(yy,2)
            for k =1:size(yy, 3)
            if yy(i,j,k) == 0 
                nanyy(i,j,k) = nan;
            end
            end
        end 
    end 
end 

function [meanMaps_, Ncc2d_, kld2d_, hi2d_, emd2d_ ]= getAccumArrays(to2dMatrix)
files = {'F:\temp\mspi\matfiles\macroMSI\24-Nov-2020_Consistent_with_Veruca_Vulgaris.mat', ...
    'F:\temp\mspi\matfiles\macroMSI\23-Nov-2020_Spitz_Nevus.mat', ...
    'F:\temp\mspi\matfiles\macroMSI\23-Nov-2020_Melanocytic_Nevus.mat', ...
    'F:\temp\mspi\matfiles\macroMSI\24-Nov-2020_Basal_Cell_Carcinoma.mat', ...
    'F:\temp\mspi\matfiles\macroMSI\24-Nov-2020_Pigmented_Seborrheic_Keratosis.mat'};

n = numel(files); 

hbMapsUnfixed = zeros(3, 3, n);
hbMapsFixed = zeros(3, 3, n);
melMapsUnfixed =zeros(3, 3, n);
melMapsFixed = zeros(3, 3, n);

hbNcc = zeros(3, 3, n);
hbKld = zeros(3, 3, n);
hbHi = zeros(3, 3, n);
hbEmd = zeros(3, 3, n);
melNcc = zeros(3, 3, n);
melKld = zeros(3, 3, n);
melHi = zeros(3, 3, n);
melEmd = zeros(3, 3, n);

for i = 1:n
    m = matfile(files{i}, 'Writable' , false);
    resultUnfixed = m.resultUnfixed;
    resultFixed = m.resultFixed; 
    roiMetricsHb = m.roiMetricsHb; 
    roiMetricsMel = m.roiMetricsMel; 
    
    hbMapsUnfixed(:,:,i) = getMeans(resultUnfixed.RoiHbMaps); 
    hbMapsFixed(:,:,i) = getMeans(resultFixed.RoiHbMaps); 
    melMapsUnfixed(:,:,i) = getMeans(resultUnfixed.RoiMelMaps);
    melMapsFixed(:,:,i) = getMeans(resultFixed.RoiMelMaps);
    
    hbNcc(:,:,i) = squeeze(roiMetricsHb(:, :, 2));
    hbKld(:,:,i) = squeeze(roiMetricsHb(:, :, 4));
    hbHi(:,:,i) = squeeze(roiMetricsHb(:, :, 3));
    hbEmd(:,:,i) = squeeze(roiMetricsHb(:, :, 5));
    
    melNcc(:,:,i) = squeeze(roiMetricsMel(:, :, 2));
    melKld(:,:,i) = squeeze(roiMetricsMel(:, :, 4));
    melHi(:,:,i) = squeeze(roiMetricsMel(:, :, 3));
    melEmd(:,:,i) = squeeze(roiMetricsMel(:, :, 5));
    
    [imgArray, nccArray, kldArray, hiArray, emdArray] = formatArrays(false, [melMapsUnfixed, melMapsFixed ; hbMapsUnfixed, hbMapsFixed], ...
        [], [melNcc; hbNcc], [melKld;hbKld], [melHi;hbHi], [melEmd;hbEmd]);
    
    writecell(files(i),'roiArrays.xlsx','Sheet',1+i,'Range','A1');
    writecell(imgArray ,'roiArrays.xlsx','Sheet',1+i,'Range','A2');
    writecell(nccArray,'roiArrays.xlsx','Sheet',1+i,'Range','A11');
    writecell(kldArray,'roiArrays.xlsx','Sheet',1+i,'Range','A20');
    writecell(hiArray,'roiArrays.xlsx','Sheet',1+i,'Range','A30');
    writecell(emdArray,'roiArrays.xlsx','Sheet',1+i,'Range','A40');


end 


hbMapsUnfixed2d = to2dMatrix(hbMapsUnfixed);
hbMapsFixed2d = to2dMatrix(hbMapsFixed);
melMapsUnfixed2d = to2dMatrix(melMapsUnfixed);
melMapsFixed2d = to2dMatrix(melMapsFixed);

hbNcc2d = to2dMatrix(hbNcc);
hbKld2d = to2dMatrix(hbKld);
hbHi2d = to2dMatrix(hbHi);
hbEmd2d = to2dMatrix(hbEmd);

melNcc2d = to2dMatrix(melNcc);
melKld2d = to2dMatrix(melKld);
melHi2d = to2dMatrix(melHi);
melEmd2d = to2dMatrix(melEmd);

hbMaps2d_ = [hbMapsUnfixed2d , hbMapsFixed2d];
melMaps2d_ = [melMapsUnfixed2d , melMapsFixed2d];
meanMaps_ = [melMaps2d_; hbMaps2d_];

Ncc2d_  = [melNcc2d; hbNcc2d ];
kld2d_  = [melKld2d; hbKld2d ];
hi2d_  = [melHi2d; hbHi2d ];
emd2d_  = [melEmd2d; hbEmd2d ];

end 

function [imgArray, nccArray, kldArray, hiArray, emdArray] = formatArrays(hasVariance, meanMaps, stdMaps, Ncc2d, kld2d, hi2d, emd2d, Ncc2dS, kld2dS, hi2dS, emd2dS)

    shortMapMethods = {'AS', 'DA', 'VI'};
    roiNames = {'HHb ROI', 'Norm ROI', 'HM ROI'}; %{Hb Norm Mel};

    imgArray = cell(size(meanMaps,1)+2, size(meanMaps,2)+2);
    imgArray{1, 3} = 'Unfixed';
    imgArray{1, 6} = 'Fixed';
    imgArray(2,3:5) = roiNames;
    imgArray(2,6:8) = roiNames;
    for i = 1:size(meanMaps,1)
        i_ = i + 2;
        if i == 1
            imgArray{i_,1} = 'Mel CrMap';
        end 
        if i == 4 
            imgArray{i_,1} = 'HbT CrMap';
        end 
        imgArray{i_,2} = shortMapMethods{mod(i-1, 3) +1};
        for j = 1:size(meanMaps,2)
            if hasVariance 
                imgArray{i_,j+2} = strcat(num2str(meanMaps(i,j), '%.3f'), char(177), num2str(stdMaps(i,j), '%.3f') );
            else 
                imgArray{i_,j+2} = num2str(meanMaps(i,j), '%.3f');
            end 
        end 
    end 

    nccArray = cell(size(Ncc2d,1)+1, size(Ncc2d,2)+2);
    kldArray = cell(size(Ncc2d,1)+1, size(Ncc2d,2)+2);
    hiArray = cell(size(Ncc2d,1)+1, size(Ncc2d,2)+2);
    emdArray = cell(size(Ncc2d,1)+1, size(Ncc2d,2)+2);

    nccArray{1, 2} = 'NCC';
    kldArray{1, 2} = 'KLD';
    hiArray{1, 2} = 'HI';
    emdArray{1, 2} = 'EMD';

    nccArray(1,3:5) = roiNames;
    kldArray(1,3:5) = roiNames;
    hiArray(1,3:5) = roiNames;
    emdArray(1,3:5) = roiNames;

    for i = 1:size(Ncc2d,1)
        i_ = i + 1; 
        if i == 1
            nccArray{i_,1} = 'Mel CrMap';
            kldArray{i_,1} = 'Mel CrMap';
            hiArray{i_,1} = 'Mel CrMap';
            emdArray{i_,1} = 'Mel CrMap';
        end 
        if i == 4 
            nccArray{i_,1} = 'HbT CrMap';
            kldArray{i_,1} = 'HbT CrMap';
            hiArray{i_,1} = 'HbT CrMap';
            emdArray{i_,1} = 'HbT CrMap';
        end 
        nccArray{i_,2} = shortMapMethods{mod(i-1, 3) +1};
        kldArray{i_,2} = shortMapMethods{mod(i-1, 3) +1};
        hiArray{i_,2} = shortMapMethods{mod(i-1, 3) +1};
        emdArray{i_,2} = shortMapMethods{mod(i-1, 3) +1};

        for j = 1:size(Ncc2d,2)
            if hasVariance 
                nccArray{i_,j+2} = strcat(num2str(Ncc2d(i,j), '%.3f'), char(177), num2str(Ncc2dS(i,j), '%.3f') );
                kldArray{i_,j+2} = strcat(num2str(kld2d(i,j), '%.3f'), char(177), num2str(kld2dS(i,j), '%.3f') );
                hiArray{i_,j+2} = strcat(num2str(hi2d(i,j), '%.3f'), char(177), num2str(hi2dS(i,j), '%.3f') );
                emdArray{i_,j+2} = strcat(num2str(emd2d(i,j), '%.3f'), char(177), num2str(emd2dS(i,j), '%.3f') );
            else 
                nccArray{i_,j+2} = num2str(Ncc2d(i,j), '%.3f');
                kldArray{i_,j+2} = num2str(kld2d(i,j), '%.3f');
                hiArray{i_,j+2} = num2str(hi2d(i,j), '%.3f');
                emdArray{i_,j+2} = num2str(emd2d(i,j), '%.3f');
 
            end 
        end 
    end 
end 

function [ ]= plotGraphs()
% files = {'F:\temp\mspi\matfiles\macroMSI\24-Nov-2020_Consistent_with_Veruca_Vulgaris.mat', ...
%     'F:\temp\mspi\matfiles\macroMSI\23-Nov-2020_Spitz_Nevus.mat', ...
%     'F:\temp\mspi\matfiles\macroMSI\23-Nov-2020_Melanocytic_Nevus.mat', ...
%     'F:\temp\mspi\matfiles\macroMSI\24-Nov-2020_Basal_Cell_Carcinoma.mat', ...
%     'F:\temp\mspi\matfiles\macroMSI\24-Nov-2020_Pigmented_Seborrheic_Keratosis.mat'};
% 
% n = numel(files); 
% 
% hbMapsUnfixed = zeros(3, 3, n);
% hbMapsFixed = zeros(3, 3, n);
% melMapsUnfixed =zeros(3, 3, n);
% melMapsFixed = zeros(3, 3, n);
% 
% hbNcc = zeros(3, 3, n);
% hbKld = zeros(3, 3, n);
% hbHi = zeros(3, 3, n);
% hbEmd = zeros(3, 3, n);
% melNcc = zeros(3, 3, n);
% melKld = zeros(3, 3, n);
% melHi = zeros(3, 3, n);
% melEmd = zeros(3, 3, n);
% 
% for i = 1:n
%     m = matfile(files{i}, 'Writable' , false);
%     resultUnfixed = m.resultUnfixed;
%     resultFixed = m.resultFixed; 
%     roiMetricsHb = m.roiMetricsHb; 
%     roiMetricsMel = m.roiMetricsMel; 
%     
%     hbMapsUnfixed(:,:,i) = getMeans(resultUnfixed.RoiHbMaps); 
%     hbMapsFixed(:,:,i) = getMeans(resultFixed.RoiHbMaps); 
%     melMapsUnfixed(:,:,i) = getMeans(resultUnfixed.RoiMelMaps);
%     melMapsFixed(:,:,i) = getMeans(resultFixed.RoiMelMaps);
%     
%     hbNcc(:,:,i) = squeeze(roiMetricsHb(:, :, 2));
%     hbKld(:,:,i) = squeeze(roiMetricsHb(:, :, 4));
%     hbHi(:,:,i) = squeeze(roiMetricsHb(:, :, 3));
%     hbEmd(:,:,i) = squeeze(roiMetricsHb(:, :, 5));
%     
%     melNcc(:,:,i) = squeeze(roiMetricsMel(:, :, 2));
%     melKld(:,:,i) = squeeze(roiMetricsMel(:, :, 4));
%     melHi(:,:,i) = squeeze(roiMetricsMel(:, :, 3));
%     melEmd(:,:,i) = squeeze(roiMetricsMel(:, :, 5));
% end 
% 
% save('spieInterimArrays.mat');

load('spieInterimArrays.mat');
load('roiArrays.mat');

roiMarker = ['r', 'g', 'b'];
methodMarker = ['^', 'o', 'd'];
methodNames = {'AS', 'DA', 'VI'};
roiNames = {'HHbT ROI', 'Norm ROI', 'HM ROI'};
mksz = 10;

%% Unfixed plot 
figure(1); clf;
h = zeros(6,1);
hold on; 
k = 0;
for j =1:3 
k = k + 1;
h(k) = plot(nan, nan, strcat('s', roiMarker(j)), 'MarkerSize', mksz, 'DisplayName', roiNames{j});
end 
for i=1:3
k = k + 1;
h(k) = plot(nan, nan, strcat(methodMarker(i), 'k'), 'MarkerSize', mksz, 'DisplayName', methodNames{i});
end

for i=1:3
    for j=1:3
        for k = 1:n
            if ~(k == 5 && j ==1)
            x = hbMapsUnfixed(i,j,k);
            y = melMapsUnfixed(i,j,k);
            plot(x, y, strcat(roiMarker(j), methodMarker(i)), 'MarkerSize', mksz);
            end
            if (k == 5 && j == 1)
                hbMapsUnfixed(i,j,k) = nan;
                melMapsUnfixed(i,j,k) = nan;
            end 
        end 
    end 
end

to2dMatrix = @to2dMatrixMean;
hbMapsUnfixed2d = to2dMatrix(hbMapsUnfixed);
melMapsUnfixed2d = to2dMatrix(melMapsUnfixed);

for i=1:3
    for j=1:3
        
        x = hbMapsUnfixed2d(i,j);
        
        y = melMapsUnfixed2d(i,j);
        plot(x, y, strcat(roiMarker(j), methodMarker(i)), 'MarkerSize', mksz, 'MarkerFaceColor', roiMarker(j));
    end 
end 
        
hold off;
xlim([0, 1]);
ylim([0, 1]);
ylabel('Mel CrMap relative concentration', 'FontSize', 13);
xlabel('HbT CrMap relative concentration', 'FontSize', 13);
legend(h, 'Location', 'eastoutside', 'FontSize', 12);
title('Average values of CrMap ROI for unfixed tissue', 'FontSize', 13);
%set(gcf,  'units','normalized','outerposition',[0 0 1 1]);
setSetting('plotName', mkNewDir('..\..\..\output\spieMaps\all\meanValues-unfixed'));
setSetting('saveImages', true);
savePlot(1);

%% Fixed plot 
figure(2); clf;
h = zeros(6,1);
hold on; 
k = 0;
for j =1:3 
k = k + 1;
h(k) = plot(nan, nan, strcat('s', roiMarker(j)), 'MarkerSize', mksz, 'DisplayName', roiNames{j});
end 
for i=1:3
k = k + 1;
h(k) = plot(nan, nan, strcat(methodMarker(i), 'k'), 'MarkerSize', mksz, 'DisplayName', methodNames{i});
end

for i=1:3
    for j=1:3
        for k = 1:n
            if ~(k == 5 && j ==1)
            x = hbMapsFixed(i,j,k);
            y = melMapsFixed(i,j,k);
            plot(x, y, strcat(roiMarker(j), methodMarker(i)), 'MarkerSize', mksz);
            end
            if (k == 5 && j == 1)
                hbMapsFixed(i,j,k) = nan;
                melMapsFixed(i,j,k) = nan;
            end
        end 
    end 
end

hbMapsFixed2d = to2dMatrix(hbMapsFixed);
melMapsFixed2d = to2dMatrix(melMapsFixed);


for i=1:3
    for j=1:3
        x = hbMapsFixed2d(i,j);
        y = melMapsFixed2d(i,j);
        plot(x, y, strcat(roiMarker(j), methodMarker(i)), 'MarkerSize', mksz, 'MarkerFaceColor', roiMarker(j));
    end 
end 
        
hold off;
xlim([0, 1]);
ylim([0, 1]);
ylabel('Mel CrMap relative concentration', 'FontSize', 13);
xlabel('HbT CrMap relative concentration', 'FontSize', 13);
legend(h, 'Location', 'eastoutside', 'FontSize', 12);
title('Average values of CrMap ROI for fixed tissue', 'FontSize', 13);
%set(gcf,  'units','normalized','outerposition',[0 0 1 1]);
setSetting('plotName', mkNewDir('..\..\..\output\spieMaps\all\meanValues-fixed'));
setSetting('saveImages', true);
savePlot(2);

    %% NCC 
    figure(1);clf;
    h = zeros(6,1);
    hold on; 
    k = 0;
    for j =1:3 
    k = k + 1;
    h(k) = plot(nan, nan, strcat('s', roiMarker(j)), 'MarkerSize', mksz, 'DisplayName', roiNames{j});
    end 
    for i=1:3
    k = k + 1;
    h(k) = plot(nan, nan, strcat(methodMarker(i), 'k'), 'MarkerSize', mksz, 'DisplayName', methodNames{i});
    end
    
    for i=1:3
        for j=1:3
            for k=1:n
                if ~(k==5 && j==1)
                    x = hbNcc(i,j,k);
                    y = melNcc(i,j,k);
                    plot(x, y, strcat(roiMarker(j), methodMarker(i)), 'MarkerSize', mksz);
                else 
                    hbNcc(i,j,k) = nan;
                    melNcc(i,j,k) = nan;
                end 
            end
        end 
    end 
    
    hbNcc2d = to2dMatrix(hbNcc);
    melNcc2d = to2dMatrix(melNcc);

    for i=1:3
        for j=1:3
            x = hbNcc2d(i,j);
            y = melNcc2d(i,j);
            plot(x, y, strcat(roiMarker(j), methodMarker(i)), 'MarkerSize', mksz,  'MarkerFaceColor', roiMarker(j));
        end 
    end 
    hold off;
    xlim([0,0.6]);
    ylim([0,0.6]);
    ylabel('NCC of Mel CrMap', 'FontSize', 13);
    xlabel('NCC of HbT CrMap', 'FontSize', 13);
    title('Normalized Correlation Coefficients', 'FontSize', 13);
    legend(h, 'Location', 'eastoutside', 'FontSize', 12);
    %set(gcf,  'units','normalized','outerposition',[0 0 1 1]);
    setSetting('plotName', '..\..\..\output\spieMaps\all\nccValues');
    setSetting('saveImages', true);
    savePlot(1);
    


    %% HI 
    figure(2);clf;
    h = zeros(6,1);
    hold on; 
    k = 0;
    for j =1:3 
    k = k + 1;
    h(k) = plot(nan, nan, strcat('s', roiMarker(j)), 'MarkerSize', mksz, 'DisplayName', roiNames{j});
    end 
    for i=1:3
    k = k + 1;
    h(k) = plot(nan, nan, strcat(methodMarker(i), 'k'), 'MarkerSize', mksz, 'DisplayName', methodNames{i});
    end
    
    for i=1:3
    for j=1:3
        for k=1:n
            if ~(k==5 && j==1)
                x = hbHi(i,j,k);
                y = melHi(i,j,k);
                plot(x, y, strcat(roiMarker(j), methodMarker(i)), 'MarkerSize', mksz);
            else 
                hbHi(i,j,k) = nan;
                melHi(i,j,k) = nan;
            end 
        end
    end 
    end 
    
    hbHi2d = to2dMatrix(hbHi);
    melHi2d = to2dMatrix(melHi);

    for i=1:3
        for j=1:3
            x = hbHi2d(i,j);
            y = melHi2d(i,j);
            plot(x, y, strcat(roiMarker(j), methodMarker(i)), 'MarkerSize', mksz,  'MarkerFaceColor', roiMarker(j));
        end 
    end 
    x = 0:0.1:1;
    y = 0:0.1:1;
    plot(x,y, 'k--');
    hold off;
    xlim([0,0.8]);
    ylim([0,0.8]);
    ylabel('HI of Mel CrMap', 'FontSize', 13);
    xlabel('HI of HbT CrMap', 'FontSize', 13);
    legend(h, 'Location', 'eastoutside', 'FontSize', 12);
    title('Histogram Intersection', 'FontSize', 13);
    %set(gcf,  'units','normalized','outerposition',[0 0 1 1]);
    setSetting('plotName', '..\..\..\output\spieMaps\all\hiValues');
    setSetting('saveImages', true);
    savePlot(2);
    
    %% EMD
    figure(3);clf;
    h = zeros(6,1);
    hold on; 
    k = 0;
    for j =1:3 
    k = k + 1;
    h(k) = plot(nan, nan, strcat('s', roiMarker(j)), 'MarkerSize', mksz, 'DisplayName', roiNames{j});
    end 
    for i=1:3
    k = k + 1;
    h(k) = plot(nan, nan, strcat(methodMarker(i), 'k'), 'MarkerSize', mksz, 'DisplayName', methodNames{i});
    end
    
    for i=1:3
    for j=1:3
        for k=1:n
            if ~(k==5 && j==1)
                x = hbEmd(i,j,k);
                y = melEmd(i,j,k);
                plot(x, y, strcat(roiMarker(j), methodMarker(i)), 'MarkerSize', mksz);
            else 
                hbEmd(i,j,k) = nan;
                melEmd(i,j,k) = nan;
            end 
        end
    end 
    end 
    
    hbEmd2d = to2dMatrix(hbEmd);
    melEmd2d = to2dMatrix(melEmd);

    for i=1:3
        for j=1:3
            x = hbEmd2d(i,j);
            y = melEmd2d(i,j);
            plot(x, y, strcat(roiMarker(j), methodMarker(i)), 'MarkerSize', mksz,  'MarkerFaceColor', roiMarker(j));
        end 
    end 
    x = 0:0.1:1;
    y = 0:0.1:1;
    plot(x,y, 'k--');
    hold off;
    xlim([0,0.5]);
    ylim([0,0.5]);
    ylabel('EMD of Mel CrMap', 'FontSize', 13);
    xlabel('EMD of HbT CrMap', 'FontSize', 13);
    legend(h, 'Location', 'eastoutside', 'FontSize', 12);
    title('Earth Mover''s Distance', 'FontSize', 13);
    %set(gcf,  'units','normalized','outerposition',[0 0 1 1]);
    setSetting('plotName', '..\..\..\output\spieMaps\all\emdValues');
    setSetting('saveImages', true);
    savePlot(3);

    %% EMD-2
    figure(4);clf;
    h = zeros(3,1);
    hold on; 
    k = 0;
    for j =1:3 
    k = k + 1;
    h(k) = plot(nan, nan, strcat('s', roiMarker(j)), 'MarkerSize', mksz, 'DisplayName', roiNames{j});
    end 
    
    for i=3
    for j=1:3
        for k=1:n
            if ~(k==5 && j==1)
                x = hbEmd(i,j,k);
                y = melEmd(i,j,k);
                plot(x, y, strcat(roiMarker(j), methodMarker(i)), 'MarkerSize', mksz);
            else 
                hbEmd(i,j,k) = nan;
                melEmd(i,j,k) = nan;
            end 
        end
    end 
    end 
    
    hbEmd2d = to2dMatrix(hbEmd);
    melEmd2d = to2dMatrix(melEmd);

    for i=3
        for j=1:3
            x = hbEmd2d(i,j);
            y = melEmd2d(i,j);
            plot(x, y, strcat(roiMarker(j), methodMarker(i)), 'MarkerSize', mksz,  'MarkerFaceColor', roiMarker(j));
        end 
    end 
    x = 0:0.1:1;
    y = 0:0.1:1;
    plot(x,y, 'k--');
    hold off;
    xlim([0,0.5]);
    ylim([0,0.5]);
    ylabel('EMD of Mel CrMap', 'FontSize', 13);
    xlabel('EMD of HbT CrMap', 'FontSize', 13);
    legend(h, 'Location', 'eastoutside', 'FontSize', 12);
    title('Earth Mover''s Distance for VI', 'FontSize', 13);
    %set(gcf,  'units','normalized','outerposition',[0 0 1 1]);
    setSetting('plotName', '..\..\..\output\spieMaps\all\emdValues-vi');
    setSetting('saveImages', true);
    savePlot(4);
    
end 

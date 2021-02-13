% for i = 1:3 
%     for j = 1:3
%     I1 = resultFixed.RoiHbMaps{i,j};
%     I2 = resultUnfixed.RoiHbMaps{i,j};
%     v1 = histogram(I1(:), 'BinEdges', linspace(0,1,21)).Values
%     v2 = histogram(I2(:), 'BinEdges', linspace(0,1,21)).Values
%     
%     %[h,p] = ttest(I1(:),I2(:),'Alpha',0.05);
%     [h,p] = ttest(v1,v2,'Alpha',0.05);
%     
%     fprintf('Rejects hypothesis: %d, P-value: %.3f.\n', h, p);
%     end 
% end 

makeExcel = false; 
if (makeExcel)
    shortMapMethods = {'AS', 'DA', 'VI'};
    roiss = {'Hb', 'Norm', 'Mel'};
    % k = 0;
    for i = 1:3 
        k = k + 1;

        Sample{k} = lesion;
        Method{k} = shortMapMethods{i};


        HMAvgUnfixedHb(k) = mean(resultUnfixed.RoiHbMaps{i,3}(:));
        HMAvgFixedHb(k) = mean(resultFixed.RoiHbMaps{i,3}(:));
        HHbAvgUnfixedHb(k) = mean(resultUnfixed.RoiHbMaps{i,1}(:));
        HHbAvgFixedHb(k) = mean(resultFixed.RoiHbMaps{i,1}(:));
        NormAvgUnfixedHb(k) = mean(resultUnfixed.RoiHbMaps{i,2}(:));
        NormAvgFixedHb(k) = mean(resultFixed.RoiHbMaps{i,2}(:));

        HMAvgUnfixedMel(k) = mean(resultUnfixed.RoiMelMaps{i,3}(:));
        HMAvgFixedMel(k) = mean(resultFixed.RoiMelMaps{i,3}(:));
        HHbAvgUnfixedMel(k) = mean(resultUnfixed.RoiMelMaps{i,1}(:));
        HHbAvgFixedMel(k) = mean(resultFixed.RoiMelMaps{i,1}(:));
        NormAvgUnfixedMel(k) = mean(resultUnfixed.RoiMelMaps{i,2}(:));
        NormAvgFixedMel(k) = mean(resultFixed.RoiMelMaps{i,2}(:));

    %     HHbStdUnfixedHb(k) = std(resultUnfixed.RoiHbMaps{i,1}(:));
    %     HHbStdFixedHb(k) = std(resultFixed.RoiHbMaps{i,2}(:));
    %     HMStdUnfixedHb(k) = std(resultUnfixed.RoiHbMaps{i,2}(:));
    %     HMStdFixedHb(k) = std(resultFixed.RoiHbMaps{i,2}(:));
    %     NormStdUnfixedHb(k) = std(resultUnfixed.RoiHbMaps{i,3}(:));
    %     NormStdFixedHb(k) = std(resultFixed.RoiHbMaps{i,3}(:));
    % 
    %     HMStdUnfixedMel(k) = std(resultUnfixed.RoiMelMaps{i,2}(:));
    %     HMStdFixedMel(k) = std(resultFixed.RoiMelMaps{i,2}(:));
    %     HHbStdUnfixedMel(k) = std(resultUnfixed.RoiMelMaps{i,1}(:));
    %     HHbStdFixedMel(k) = std(resultFixed.RoiMelMaps{i,2}(:));
    %     NormStdUnfixedMel(k) = std(resultUnfixed.RoiMelMaps{i,3}(:));
    %     NormStdFixedMel(k) = std(resultFixed.RoiMelMaps{i,3}(:));


        HMRoiMetricsHbNCC(k) = roiMetricsHbTable.NCC(7 + i - 1);
        HMRoiMetricsMelNCC(k) = roiMetricsMelTable.NCC(7 + i - 1);
        HHbRoiMetricsHbNCC(k) = roiMetricsHbTable.NCC(1 + i - 1);
        HHbRoiMetricsMelNCC(k) = roiMetricsMelTable.NCC(1 + i - 1);
        NormRoiMetricsHbNCC(k) = roiMetricsHbTable.NCC(4 + i - 1);
        NormRoiMetricsMelNCC(k) = roiMetricsMelTable.NCC(4 + i - 1);
        MetricsHbNCC(k) = metricsHbTable.NCC(i);
        MetricsMelNCC(k) = metricsMelTable.NCC(i);

        HMRoiMetricsHbHI(k) = roiMetricsHbTable.HI(7 + i - 1);
        HMRoiMetricsMelHI(k) =  roiMetricsMelTable.HI(7 + i - 1);
        HHbRoiMetricsHbHI(k) = roiMetricsHbTable.HI(1 + i - 1);
        HHbRoiMetricsMelHI(k) = roiMetricsMelTable.HI(1 + i - 1);
        NormRoiMetricsHbHI(k) = roiMetricsHbTable.HI(4 + i - 1);
        NormRoiMetricsMelHI(k) = roiMetricsMelTable.HI(4 + i - 1);
        MetricsHbHI(k) = metricsHbTable.HI(i);
        MetricsMelHI(k) = metricsMelTable.HI(i);

        HMRoiMetricsHbKLD(k) =  roiMetricsHbTable.KLD(7 + i - 1);
        HMRoiMetricsMelKLD(k) =  roiMetricsHbTable.KLD(7 + i - 1);
        HHbRoiMetricsHbKLD(k) =  roiMetricsHbTable.KLD(1 + i - 1);
        HHbRoiMetricsMelKLD(k) = roiMetricsHbTable.KLD(1 + i - 1);
        NormRoiMetricsHbKLD(k) =  roiMetricsHbTable.KLD(4 + i - 1);
        NormRoiMetricsMelKLD(k) =  roiMetricsHbTable.KLD(4 + i - 1);
        MetricsHbKLD(k) =   metricsHbTable.KLD(i);
        MetricsMelKLD(k) =  metricsMelTable.KLD(i);
    end 
end 

%% Make images for mean values 
makeImage = true; 
oneFigure = false;

if makeImage
    close all; 
    load('roiArrays.mat');

    methodMarker = ['r', 'g', 'b'];
    roiMarker = ['^', 'o', 'd'];
    methodNames = {'AS', 'DA', 'VI'};
    roiNames = {'HHbT ROI', 'Norm ROI', 'HM ROI'};
    mksz = 10;
        
    if oneFigure
        %% for single figure 
        
        figure(1); clf;
        hold on; 
        k = 0;
        h = zeros(8,1);
        for j =1:3 
        k = k + 1;
        h(k) = plot(nan, nan, strcat('s', methodMarker(j)), 'MarkerSize', mksz, 'DisplayName', methodNames{j});
        end 
        for i=1:3
        k = k + 1;
        h(k) = plot(nan, nan, strcat(roiMarker(i), 'k'), 'MarkerSize', mksz, 'DisplayName', roiNames{i});
        end
        k =k+1;
        h(k) = plot(nan,nan, strcat('s', 'k'), 'MarkerSize', mksz, 'DisplayName', 'Unfixed tissue');
        k = k+1;
        h(k) = plot(nan,nan, strcat('s', 'k'), 'MarkerSize', mksz, 'DisplayName', 'Fixed tissue', 'MarkerFaceColor', 'k');
        
        for i=1:3
            for j=1:3
                tmp = strsplit(imgArray{5+i,2+j}, 'Å}');
                x = str2double(tmp{1});
                tmp = strsplit(imgArray{2+i,2+j}, 'Å}');
                y = str2double(tmp{1});
                plot(x, y, strcat(roiMarker(j), methodMarker(i)), 'MarkerSize', mksz);
            end 
            for j=1:3
                tmp = strsplit(imgArray{5+i,5+j}, 'Å}');
                x = str2double(tmp{1});
                tmp = strsplit(imgArray{2+i,5+j}, 'Å}');
                y = str2double(tmp{1});
                plot(x, y, strcat(roiMarker(j), methodMarker(i)), 'MarkerSize', mksz, 'MarkerFaceColor', methodMarker(i));
            end
        end 
        hold off;
        xlim([0, 0.6]);
        ylim([0, 1]);
        ylabel('Mel CrMap relative concentration', 'FontSize', 10);
        xlabel('HbT CrMap relative concentration', 'FontSize', 10);
        legend(h, 'Location', 'eastoutside', 'FontSize', 12);
        %set(gcf,  'units','normalized','outerposition',[0 0 1 1]);
        setSetting('plotName', '..\..\..\output\spieMaps\meanValues');
        setSetting('saveImages', true);
        savePlot(1);
    end 
    
    if ~oneFigure
        %% for unfixed 
        figure(1); clf;
        h = zeros(6,1);
        hold on; 
        k = 0;
        for j =1:3 
        k = k + 1;
        h(k) = plot(nan, nan, strcat('s', methodMarker(j)), 'MarkerSize', mksz, 'DisplayName', methodNames{j});
        end 
        for i=1:3
        k = k + 1;
        h(k) = plot(nan, nan, strcat(roiMarker(i), 'k'), 'MarkerSize', mksz, 'DisplayName', roiNames{i});
        end
                
        for i=1:3
            for j=1:3
                tmp = strsplit(imgArray{5+i,2+j}, 'Å}');
                x = str2double(tmp{1});
                tmp = strsplit(imgArray{2+i,2+j}, 'Å}');
                y = str2double(tmp{1});
                plot(x, y, strcat(roiMarker(j), methodMarker(i)), 'MarkerSize', mksz);
            end 
        end 
        hold off;
        xlim([0, 0.6]);
        ylim([0, 1]);
        ylabel('Mel CrMap relative concentration', 'FontSize', 10);
        xlabel('HbT CrMap relative concentration', 'FontSize', 10);
        legend(h, 'Location', 'eastoutside', 'FontSize', 12);
        %set(gcf,  'units','normalized','outerposition',[0 0 1 1]);
        setSetting('plotName', '..\..\..\output\spieMaps\meanValues-unfixed');
        setSetting('saveImages', true);
        savePlot(1);
        
        %% for fixed 
        figure(2); clf;
        hold on; 
        h = zeros(6,1);
        k = 0;
        for j =1:3 
        k = k + 1;
        h(k) = plot(nan,nan, strcat('s', methodMarker(j)), 'MarkerSize', mksz, 'DisplayName', methodNames{j}, 'MarkerFaceColor', methodMarker(j));
        end 
        for i=1:3
        k = k + 1;
        h(k) = plot(nan,nan, strcat(roiMarker(i), 'k'), 'MarkerSize', mksz, 'DisplayName', roiNames{i}, 'MarkerFaceColor', 'k');
        end
        
        for i=1:3
            for j=1:3
                tmp = strsplit(imgArray{5+i,5+j}, 'Å}');
                x = str2double(tmp{1});
                tmp = strsplit(imgArray{2+i,5+j}, 'Å}');
                y = str2double(tmp{1});
                plot(x, y, strcat(roiMarker(j), methodMarker(i)), 'MarkerSize', mksz, 'MarkerFaceColor', methodMarker(i));
            end 
        end 
        hold off;
        xlim([0,0.6]);
        ylim([0,1]);
        ylabel('Mel CrMap relative concentration', 'FontSize', 10);
        xlabel('HbT CrMap relative concentration', 'FontSize', 10);
        legend(h, 'Location', 'eastoutside', 'FontSize', 12);
        %set(gcf,  'units','normalized','outerposition',[0 0 1 1]);
        setSetting('plotName', '..\..\..\output\spieMaps\meanValues-fixed');
        setSetting('saveImages', true);
        savePlot(2);
        
    end 
end 


%% Make images for metrics 
makeImages2 = false; 
if makeImages2 
    close all; 
    load('roiArrays.mat');

    methodMarker = ['r', 'g', 'b'];
    roiMarker = ['^', 'o', 'd'];
    methodNames = {'AS', 'DA', 'VI'};
    roiNames = {'HHbT ROI', 'Norm ROI', 'HM ROI'};
    mksz = 10;
    
    %% NCC 
    figure(1);clf;
    hold on; 
    h = zeros(6,1);
    k = 0;
    for j =1:3 
    k = k + 1;
    h(k) = plot(nan,nan, strcat('s', methodMarker(j)), 'MarkerSize', mksz, 'DisplayName', methodNames{j}, 'MarkerFaceColor', methodMarker(j));
    end 
    for i=1:3
    k = k + 1;
    h(k) = plot(nan,nan, strcat(roiMarker(i), 'k'), 'MarkerSize', mksz, 'DisplayName', roiNames{i}, 'MarkerFaceColor', 'k');
    end
    
    for i=1:3
        for j=1:3
            tmp = strsplit(nccArray{4+i,2+j}, 'Å}');
            x = str2double(tmp{1});
            tmp = strsplit(nccArray{1+i,2+j}, 'Å}');
            y = str2double(tmp{1});
            plot(x, y, strcat(roiMarker(j), methodMarker(i)), 'MarkerSize', mksz);
        end 
    end 
    hold off;
    xlim([0,0.4]);
    ylim([0,0.4]);
    ylabel('NCC of Mel CrMap', 'FontSize', 10);
    xlabel('NCC of HbT CrMap', 'FontSize', 10);
    legend(h, 'Location', 'eastoutside', 'FontSize', 12);
    %set(gcf,  'units','normalized','outerposition',[0 0 1 1]);
    setSetting('plotName', '..\..\..\output\spieMaps\nccValues');
    setSetting('saveImages', true);
    savePlot(1);
        
    %% HI 
    figure(2);clf;
    hold on; 
    h = zeros(6,1);
    k = 0;
    for j =1:3 
    k = k + 1;
    h(k) = plot(nan,nan, strcat('s', methodMarker(j)), 'MarkerSize', mksz, 'DisplayName', methodNames{j}, 'MarkerFaceColor', methodMarker(j));
    end 
    for i=1:3
    k = k + 1;
    h(k) = plot(nan,nan, strcat(roiMarker(i), 'k'), 'MarkerSize', mksz, 'DisplayName', roiNames{i}, 'MarkerFaceColor', 'k');
    end
    
    for i=1:3
        for j=1:3
            tmp = strsplit(hiArray{4+i,2+j}, 'Å}');
            x = str2double(tmp{1});
            tmp = strsplit(hiArray{1+i,2+j}, 'Å}');
            y = str2double(tmp{1});
            plot(x, y, strcat(roiMarker(j), methodMarker(i)), 'MarkerSize', mksz);
        end 
    end 
    hold off;
    xlim([0,0.8]);
    ylim([0,0.8]);
    ylabel('HI of Mel CrMap', 'FontSize', 10);
    xlabel('HI of HbT CrMap', 'FontSize', 10);
    legend(h, 'Location', 'eastoutside', 'FontSize', 12);
    %set(gcf,  'units','normalized','outerposition',[0 0 1 1]);
    setSetting('plotName', '..\..\..\output\spieMaps\hiValues');
    setSetting('saveImages', true);
    savePlot(2);
    
    %% EMD
    figure(3);clf;
    hold on; 
    h = zeros(6,1);
    k = 0;
    for j =1:3 
    k = k + 1;
    h(k) = plot(nan,nan, strcat('s', methodMarker(j)), 'MarkerSize', mksz, 'DisplayName', methodNames{j}, 'MarkerFaceColor', methodMarker(j));
    end 
    for i=1:3
    k = k + 1;
    h(k) = plot(nan,nan, strcat(roiMarker(i), 'k'), 'MarkerSize', mksz, 'DisplayName', roiNames{i}, 'MarkerFaceColor', 'k');
    end
    
    for i=1:3
        for j=1:3
            tmp = strsplit(emdArray{4+i,2+j}, 'Å}');
            x = str2double(tmp{1});
            tmp = strsplit(emdArray{1+i,2+j}, 'Å}');
            y = str2double(tmp{1});
            plot(x, y, strcat(roiMarker(j), methodMarker(i)), 'MarkerSize', mksz);
        end 
    end 
    hold off;
    xlim([0,0.4]);
    ylim([0,0.4]);
    ylabel('EMD of Mel CrMap', 'FontSize', 10);
    xlabel('EMD of HbT CrMap', 'FontSize', 10);
    legend(h, 'Location', 'eastoutside', 'FontSize', 12);
    %set(gcf,  'units','normalized','outerposition',[0 0 1 1]);
    setSetting('plotName', '..\..\..\output\spieMaps\emdValues');
    setSetting('saveImages', true);
    savePlot(3);
    
end 


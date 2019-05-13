function [] = plots(plotType, varargin)
%% PLOTS Each time plots an informative plot based on 'plotType'
%
% Syntax:
%   plots(plotType, 1, curves, 'leftcamera23')
% Additional (Name,Value) arguments are available

%%%%%%%%%%%%%%%%%%%%%%%%Initialize%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
defaultfc = [450, 465, 505, 525, 575, 605, 630];
addRequired(p, 'plotType', @(x) any(validatestring(x, {'sensitivity', 'illuminationAndSensitivity', 'illumination', 'estimationComparison', ...
    'allEstimations', 'singlemeasurement', 'd65', 'methodErrors', 'overlapSpectrum', 'overlapSpectrumSample', 'pca', 'lda', 'pca b', 'lda b', ...
    'pcalda', 'classificationErrors', 'cropped', 'segmentation', 'roc', 'performanceComparison', 'visual','classificationPerformance', 'lbp', ...
    'normSensitivity', 'featureImportance', 'classificationPerformanceBars', 'reconstructionPerformanceBars'})));
addOptional(p, 'fig', -1);
addOptional(p, 'curves', []);
addOptional(p, 'title', [],  @(x) ischar(x) || iscell(x));
addParameter(p, 'Wavelength', []);
addParameter(p, 'Sensitivity', []);
addParameter(p, 'Illumination', []);
addParameter(p, 'Fc', defaultfc);
addParameter(p, 'PlotName', [], @(x) ischar(x));
addParameter(p, 'LineNames', {});
addParameter(p, 'Markers', {});
addParameter(p, 'Errors', []);
addParameter(p, 'Latent', []);
addParameter(p, 'Explained', []);
addParameter(p, 'SaveOptions', []);
addParameter(p, 'Image', []);
addParameter(p, 'Coordinates', []);
addParameter(p, 'Performance', []);
addParameter(p, 'FoldPerformance', []);
addParameter(p, 'Cmap', 'jet');
addParameter(p, 'Overlay', []);
addParameter(p, 'AdditionalImage', []);
addParameter(p, 'Alpha', 0.5);
addParameter(p, 'saveInBW', false);

parse(p, plotType, varargin{:});
plotType = p.Results.plotType;
fig = p.Results.fig;
curves = p.Results.curves;
figTitle = p.Results.title;
sensitivity = p.Results.Sensitivity;
wavelength = p.Results.Wavelength;
illumination = p.Results.Illumination;
fc = p.Results.Fc;
plotName = p.Results.PlotName;
lineNames = p.Results.LineNames;
markers = p.Results.Markers;
errors = p.Results.Errors;
latent = p.Results.Latent;
explained = p.Results.Explained;
saveOptions = p.Results.SaveOptions;
I = p.Results.Image;
coordinates = p.Results.Coordinates;
performance = p.Results.Performance;
foldPerformance= p.Results.FoldPerformance;
cmap = p.Results.Cmap;
overlay = p.Results.Overlay;
I2 = p.Results.AdditionalImage;
alpha = p.Results.Alpha;
saveInBW = p.Results.saveInBW;

if fig < 0
    figure;
else
    figure(fig);
    clf(fig);
end

if ~isempty(saveOptions)
    savePlot = saveOptions.saveImages;
    plotName = saveOptions.plotName;
    saveInHQ = saveOptions.saveInHQ;
else
    savePlot = false;
    saveInHQ = false;
end

if ~isempty(figTitle);  figTitle = strrep(figTitle, '_', ' '); end

w = warning('off', 'all');
switch plotType
    case {'sensitivity', 'normSensitivity', 'illumination', 'illuminationAndSensitivity'}
%% Plot camera sensitivity and illumination %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %{
            % Arguments  ( wavelength, bandWavelength, illumination )
            % wavelength: 401x1 or 81x1
            % bandWavelength: the names of the bands
            % illumination: the illumination light 401x7
        %}
        if contains(lower(plotType), 'sensitivity')
            sensitivityN = size(sensitivity, 2);
            if (saveInBW)
                color(1,:)=HInt2RGB(1,100); % red, darkest 
                color(2,:)=HInt2RGB(3,64); % green, less dark 
                color(3,:)=HInt2RGB(7,10); % blue cyan, lightest
            else 
                color = hsv(3);
                %color = flip(color);
            end

            if strcmp(plotType, 'normSensitivity'); n  = max(sensitivity(:)); else; n = 1; end
            
            hold on;
            if (sensitivityN == 7)
                fnames = cellstr(num2str(defaultfc'));
            else
                fnames = {'red'; 'green'; 'blue'};
            end
            for j = 1:sensitivityN
                plot(wavelength, sensitivity(:, j) / n, 'Color', color(j,:), 'LineWidth', 3, 'DisplayName', fnames{j});
            end
            hold off;
            if strcmp(plotType, 'normSensitivity'); ylabel('Normalized Sensitivity', 'FontSize', 17); else; ylabel('Sensitivity', 'FontSize', 17); end     
            xlabel('Wavelength \lambda (nm)', 'FontSize', 15);
        end
        
        if (contains(lower(plotType), 'sensitivity') && contains(lower(plotType), 'illumination')); yyaxis right; end

        if contains(lower(plotType), 'illumination')
            if (saveInBW)
                color(7,:)=HInt2RGB(1,100); % red, darkest 
                color(6,:)=HInt2RGB(7,82); % cyan, less dark 
                color(5,:)=HInt2RGB(3,64); % green, less dark 
                color(4,:)=HInt2RGB(5,53); %
                color(3,:)=HInt2RGB(9,46); % magenta, less dark 
                color(2,:)=HInt2RGB(2,28); % orange, less dark 
                color(1,:)=HInt2RGB(7,10); % blue cyan, lightest
            else 
                color = jet(7);
            end

            hold on;
            for j = 1:numel(fc)
                plot(wavelength, illumination(:, j) * 10, 'DisplayName', [ num2str(fc(j)), ' nm'], 'Color', color(j,:), 'LineWidth', 2);
            end
            plot(wavelength, illumination(:,8) * 10, 'DisplayName', 'white', 'LineStyle', ':', 'LineWidth', 3);
            hold off;
            xlabel('Wavelength \lambda (nm)', 'FontSize', 15);
            ylabel('Luminous Intensity', 'FontSize', 17);
        end
        
        if (contains(lower(plotType), 'sensitivity') && contains(lower(plotType), 'illumination'))
            title('Overlap of illumination and camera sensitivity spectrum')
        elseif contains(lower(plotType), 'normsensitivity') 
            title('Normalized Camera Sensitivities', 'FontSize', 15);
        elseif contains(lower(plotType), 'sensitivity') 
            title('Camera Sensitivities', 'FontSize', 15);
        elseif contains(lower(plotType), 'illumination') 
            title('Normalized Luminous Intensity', 'FontSize', 15);
        end

        ax = gca;
        ax.FontSize = 15; 
        xlim([400, 700]);  
        xticks([400, 500, 600, 700]);
        yticks([0, 0.5, 1]);
        l = legend;
        l.FontSize = 13;
        % End Plot camera sensitivity        
        
    case 'estimationComparison'       
%% Plot estimation comparison%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        %{
            % returns a comparative plot of the estimation curves versus the measured spectrum.
            %wavelength: vector, e.g. 401x1
            %curves: the estimation results in order { 'Measured', 'Est-MSgreen', 'Est-MSrms','Est-MSadjusted', 'Est-MSextended', 'Est-RGB'}
            %name: the name of the sample, is included in the plot title
            %method: the estimation method (eg SmoothingMatrixMethod)
            %plotOutName = if it's not empty, then the plot is saved with that filename
            % MSIreflectances = presents also a subplot os the MSI reflectances value
        %}
        
        % each column of 'spectrum' is that data for a plot line
        curveN = size(curves, 2);
        if isempty(lineNames)
            lineNames = {'MS center \lambda', 'Measured', 'Est-MSgreen', 'Est-MSrms', 'Est-MSadjusted', 'Est-MSextended', 'Est-RGB'};
        else
            lineNames = cellfun(@(x) strrep(x, '_', ' '), lineNames, 'un', 0);
        end
        
        if isempty(markers)
            for i = 1:curveN
                markers{i} = 'none';
            end
        end
        
        if ~exist('MSIreflectances', 'var'); plotN = 2; else; plotN = 3; end
        
        color = colorcube(curveN+10);         
        [~, peakIdx] = ismember(fc, wavelength); % mark filter wavelengths
        
        subplot(1, plotN, [1, 2]);
        hold on
        for i = 1:curveN
            if  contains( lower(lineNames{i + 1}), 'rgb')
                lineStyle = ':';
            elseif contains( lower(lineNames{i + 1}), 'measured')
                lineStyle = '--';
            else 
                lineStyle = '-';
            end
            plot(wavelength, curves(:, i) .* 100, 'Color', color(i, :), 'Marker', markers{i}, ...
                'LineWidth', 1.3, 'LineStyle', lineStyle, 'DisplayName', lineNames{i + 1}); % plot estimated reflectances
        end
        plot(wavelength(peakIdx), curves(peakIdx, 1) .* 100, 'rx', 'DisplayName', lineNames{1}, 'LineWidth', 1.3); % plot measured reflectance
        hold off
        
        xlabel('Wavelength \lambda (nm)');
        ylabel('Reflectance %');
        xlim([400, 700]);
        figTitle = strjoin({'Comparative plot of Wiener estimation for Sample', simpleSampleName(figTitle)}, ' ');
        title(figTitle);
        legend({lineNames{2:(curveN+1)}, lineNames{1}}, 'Location', 'best', 'FontSize', 12); % 'Orientation','horizontal');
        
        if (plotN > 2)
            subplot(1, plotN, 3);
            
            gg = raw2msi(MSIreflectances, 'rms');
            ref = mean(mean(gg, 3), 2);
            plot(fc, ref', 'm*-'); %magenda new
            xlim([400, 750]);
            
            title('MSI reflectances');
            xlabel('Wavelength \lambda (nm)');
            ylabel('MSI double value (rms of the 3 RGB values)')
        end
        % End Plot estimation comparison
        
    case 'allEstimations'       
%% plots all estimated curves for all pixels in an image area%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        hold on
        for i = 2:length(curves) - 1
            plot(wavelength, curves(i, :), 'm');
        end
        h1 = plot(wavelength, curves(1, :), 'b', 'LineWidth', 6); % measured
        h2 = plot(wavelength, curves(end, :), 'y', 'LineWidth', 6); % estimated
        hold off
        
        figTitle = {sprintf('Comparative plot of Wiener estimation results for a pixel area '); ...
                    sprintf('Sample %s', figTitle)}; %sprintf('Sample %s with %s smoothing matrix ', name, method)};
        figTitle = strrepAll(figTitle);
        title(figTitle);
        xlabel('Wavelength \lambda (nm)');
        ylabel('Reflectance');
        legend(h1, {'Measured reflectance'});
        legend(h2, {'Estimated reflectance'});
        % end plots all estimated curves for all pixels in an image area
             
    case 'd65'      
%% plot d65 illumination %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [lambda, illum] = illuminant('d65');
        plot(lambda, illum/max(illum), 'm', 'LineWidth', 2);
        xlabel('wavelength (nm)')
        xlim([300, 830])
        ylabel('relative spectral power distribution')
        title('The CIE D65 daylight illuminant ')
        %end plot d65 illumination
              
    case 'methodErrors'     
%% plot method errors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if contains(figTitle, 'nmse', 'IgnoreCase', true)
            avge = [errors.avgrmse];
            maxe = [errors.maxrmse];
            mine = [errors.minrmse];
            stde = [errors.stdrmse];
            labely = 'NMSE';
        else
            avge = [errors.avgnmse];
            maxe = [errors.maxnmse];
            mine = [errors.minnmse];
            stde = [errors.stdnmse];
            labely = 'RMSE';
        end
        
        if contains(figTitle, 'System') 
            pvsms = unique({errors.pixelValueSelectionMethod}, 'stable');
        end
        
        if contains(figTitle, 'Noise') 
            pvsms = unique({errors.noiseType}, 'stable');
        end
        
        if contains(figTitle, 'Matrix')
            smms = unique({errors.smoothingMatrixMethod}, 'stable');
        end 
        labelx = 'Wiener estimation variation';

        
        marker = ['o', 's', 'd', '^', '*', 'h', 'p', 'v', '<', '+', '>','o', 's', 'd', '^', '*'];
        if ((length(smms) > 1) && (length(pvsms) > 1))
            pvsmColor = parula(length(pvsms)+1);
            pvsmColor = pvsmColor(1:length(pvsms),:);
        else
            pvsmColor = jet(length(pvsms));
        end
%         pvsmColor = pvsmColor(2:9,:); % remove yellow and black color 
        
        hold on
        % Generate dummy info for plot handles "h"
        h = zeros(length(pvsms), 1);
        for i = 1:length(pvsms)
            h(i) = plot([NaN, NaN], 'Color', pvsmColor(i, :), 'Marker', marker(i), 'DisplayName', pvsms{i});
        end
        if contains(figTitle, 'avg', 'IgnoreCase', true)
            i = i + 1;
            h(i) = plot([NaN, NaN], 'ko', 'DisplayName', 'Mean and Std');
        end
        if contains(figTitle, 'max', 'IgnoreCase', true)
            i = i + 1;
            h(i) = plot([NaN, NaN], 'k-.', 'DisplayName', 'Max value');
        end
        if contains(figTitle, 'min', 'IgnoreCase', true)
            i = i + 1;
            h(i) = plot([NaN, NaN], 'k:', 'DisplayName', 'Min value');
        end
        % end dummy data for legend entries
        
        for i = 1:length(pvsms)
            for j = 1:length(smms)
                if contains(figTitle, 'avg', 'IgnoreCase', true)
                    errorbar(j, avge((i - 1)*length(smms)+j), stde((i - 1)*length(smms)+j) ./ size(stde,2) , ...
                        'Color', pvsmColor(i, :), 'LineWidth', 1, 'Marker', marker(i), 'MarkerSize', 10);
                end
            end
            if contains(figTitle, 'max', 'IgnoreCase', true)
                plot(1:length(smms), maxe, 'Color', pvsmColor(i, :), 'LineStyle', '-.');
            end
            if contains(figTitle, 'min', 'IgnoreCase', true)
                plot(1:length(smms), mine, 'Color', pvsmColor(i, :), 'LineStyle', ':');
            end
        end
        hold off
        
        grid on
        grid minor
        xlabel(labelx);
        xlim([0, numel(smms) + 1]);
        xticks(1:numel(smms));
        xticklabels(strrep(smms, '_', ' '));
        xtickangle(45);
        ylabel(labely);
        ax = gca;
        ax.YRuler.Exponent = 0;
        legend(h, 'Location', 'best');
        title('Comparison of Mean and Standard Error Values for various Estimation Configurations')
        set(gcf, 'Position', get(0, 'Screensize'));
        saveInHQ = true;
        %end  plot method errors
               
    case 'classificationErrors'       
%% Classication errors%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        n = min(20, length(errors));
        accuracy = [errors(1:n).Accuracy];
        typeI = [errors(1:n).TypeI];
        typeII = [errors(1:n).TypeII];
        hold on
        scatter(1:n, accuracy, 'rx')
        scatter(1:n, typeI, 'bo')
        scatter(1:n, typeII, 'mo')
        text(1:n, accuracy-2, arrayfun(@(x) sprintf('%.2f%%', x), accuracy, 'UniformOutput', false), 'FontSize', 8);
        text(1:n, typeI+2, arrayfun(@(x) sprintf('%.2f%%', x), typeI, 'UniformOutput', false), 'FontSize', 8);
        text(1:n, typeII-2, arrayfun(@(x) sprintf('%.2f%%', x), typeII, 'UniformOutput', false), 'FontSize', 8);
        
        hold off
        xlim([0, n + 1])
        grid on
        xticks(1:n)
        ll = cell(n, 1);
        for i = 1:n
            ll{i} = strcat(errors(i).Input, '|', errors(i).Projection, '|', num2str(errors(i).Neighbours), '-', errors(i).VoteRule, '|', errors(i).Distance);
        end
        xticklabels(ll)
        xtickangle(45)
        legend('Accuracy', 'False positive rate', 'False negative rate', 'Location', 'best')
        ylabel('Classification Accuracy')
        title(strcat('Classification results (', errors(i).Validation, ' validation)'));
                
    case {'PCA', 'LDA', 'PCA b', 'LDA b', 'PCALDA'}     
%% plot discriminant analysis results%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if strcmp(plotType, 'PCA') && ~isempty(latent) && ~isempty(explained)
            subplot(1, 3, 1);
            plot(1:length(latent), latent, '-mx');
            xlabel('Sorted eigenvalue index');
            ylabel('PCA eigenvalues');
            xlim([1, length(latent) + 1])
            text(1:length(explained), latent+10^(ceil(log10(max(latent))) - 2), arrayfun(@(x) sprintf('%.2f%%', x), explained, 'UniformOutput', false), 'FontSize', 8);
            title('Largest eigenvalues and explained covariance percentage')
            
            subplot(1, 3, [2, 3]);
        end
        
        marker = ['o', 'x', 'd', '^', '*', 'h', 'p', 'v', 's', '<', '+', '>'];
        observations = size(curves, 1);
        attr = split(lineNames, ' ');
        sample = {attr{:, 1}}';
        type = {attr{:, 2}}';
        isBenign = {attr{:, 3}}';
        idx = strcmp(isBenign, 'Normal') | strcmp(isBenign, 'Benign');
        colors(idx) = 'b';
        colors(~ismember(idx, 1:length(isBenign))) = 'r';
        
        %dummy legends
        h = zeros(2, 1);
        hold on
        h(1) = plot([NaN, NaN], 'Color', 'r', 'DisplayName', 'Malignant');
        h(2) = plot([NaN, NaN], 'Color', 'b', 'DisplayName', 'Benign');
        hold off
        
        if contains(figTitle, 'Fix', 'IgnoreCase', true)
            hold on
            h(3) = plot([NaN, NaN], 'Color', 'k', 'Marker', marker(1), 'DisplayName', 'Fixed');
            h(4) = plot([NaN, NaN], 'Color', 'k', 'Marker', marker(2), 'DisplayName', 'Unfixed');
            hold off
            markers = repmat(marker(1), observations, 1);
            idx = arrayfun(@(x) any(strcmp(x, 'unfixed')), type);
            markers(idx) = marker(2);
        end
        
        if contains(figTitle, 'Sample', 'IgnoreCase', true)
            markers = repmat(marker(1), observations, 1);
            samples = unique(sample, 'stable');
            for i = 1:length(samples)
                hold on
                h(i+2) = plot([NaN, NaN], 'Color', 'k', 'Marker', marker(i), 'DisplayName', samples{i});
                hold off
                idx = arrayfun(@(x) any(strcmp(x, samples{i})), sample);
                markers(idx) = marker(i);
            end
        end
        %dummy legends
        
        hold on
        for i = 1:observations
            scatter(curves(i, 1), curves(i, 2), [], colors(i), markers(i));
        end
        hold off
        
        if strcmp(plotType, 'pca')
            title('PCA scores for PC1&2')
            xlabel('Principal component 1');
            ylabel('Principal component 2');
        else
            title('LDA projections for LD1&2')
            xlabel('Lidear discriminant 1');
            ylabel('Linear Discriminant 2');
        end
        ax = gca;
        ax.XRuler.Exponent = 0;
        ax.YRuler.Exponent = 0;
        titl = strsplit(plotName, '\');
        figTitle = strrepAll(titl{end});
        suptitle(figTitle);
        
        legend(h, 'Location', 'best');
        set(gcf, 'Position', get(0, 'Screensize'));
        saveInHQ = true;
        
    case 'segmentation'
%% Plot Segmentation results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        h1 = subplot(1,2,1);
        clf(h1)
        imoverlay_medical(rgb2gray(I), overlay, [], [], 'hot', [], h1);
        c = colorbar('location', 'westoutside','Ticks', linspace(0,1,5), 'TickLabels', linspace(0,100,5));
        c.Label.String = 'Channel Agreement (%)';
        c.Label.FontSize = 15;
        c.Label.FontWeight = 'bold';
        c.Limits = [0,1];
        figTitle = strjoin({'Channel ROI masks for sample', strrepAll(plotName)}, ' ');
        title(figTitle, 'FontSize', 15)
        
        h2 = subplot(1,2,2);
        colors = jet(size(coordinates,1));
        imshow(I2);
        hold on
        for i = 1:size(coordinates, 1)
            x = coordinates(i,1);
            y = coordinates(i,2);
            plot(x, y, '*', 'LineWidth', 2, 'MarkerSize', 7, 'MarkerFaceColor', colors(i,:));
        end
        hold off
        figTitle = strjoin({'Final ROI of sample', strrepAll(plotName)}, ' ');
        title(figTitle, 'FontSize', 15);
        plotName = strcat(plotName, '_segments');
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 1, 1]);
                
    case 'roc'
        %% Show ROC performance of classifier %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        c = colormap(lines);    
        hold on
        folds = length(foldPerformance);
        for i = 1:folds*(folds < 11)
            if ~isempty(foldPerformance(i).AUC)
                plot(foldPerformance(i).ROCX, foldPerformance(i).ROCY, 'DisplayName', sprintf('Fold %d (AUC = %.3f)', i, foldPerformance(i).AUC), 'LineWidth', 2);
            end
        end
        plot(0:0.1:1, 0:0.1:1, 'k--', 'DisplayName', 'Chance');
        if ~isempty(performance.ROCX)
            plot(performance.ROCX(:,1), performance.ROCY(:,1), 'm-*',  'DisplayName', sprintf('Average (AUC = %.3f)', performance.AUC(1)), 'LineWidth', 2);
%         shadedErrorBar(performance.ROCX(:,1), performance.ROCY(:,1), abs(performance.ROCY(:,1) - performance.ROCY(:,2:3)),'lineprops','m-*');
        end
        hold off 
        xlim([-0.1, 1.1]);
        ylim([-0.1, 1.1]);
        xlabel('False positive rate');
        ylabel('True positive rate');
        h = findobj(gca,'Type','line');
        hh = h(contains({h.DisplayName}, {'AUC', 'Chance'}));
        legend(hh, 'Location', 'eastoutside');
        title(figTitle);
        set(gcf, 'Position', get(0, 'Screensize'));
        
    case 'LBP'
%% Plot LBP values on image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        h1 = subplot(1,2,1);
        im = imagesc(imrotate(I,90));
        colormap('hot');
        im.AlphaData = .8;
        c = colorbar;
        axis off;
        c.Label.FontSize = 20;
        c.Label.FontWeight = 'bold';
        c.Label.String = 'Normalized Texture Descriptor Value';
        title('Conventional LBP descriptor values for RGB','FontSize', 20);
        set(gcf, 'Position', get(0, 'Screensize'));
        
        h2 = subplot(1,2,2);
        im = imagesc(imrotate(I2,90));
        colormap('hot');
        im.AlphaData = .8;
        c = colorbar;
        axis off;
        c.Label.FontSize = 20;
        c.Label.FontWeight = 'bold';
        c.Label.String = 'Normalized Texture Descriptor Value';
        title('SumLBP descriptor values for MSI','FontSize', 20);
        set(gcf, 'Position', get(0, 'Screensize'));

    case 'performanceComparison'
%% Performace comparison for classifiers%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
        N = length(lineNames);
        x = 1:N;
        auc = performance(1,:)';
        accur = performance(2,:)';
        hold on 

        if ~isnan(auc)
            yyaxis left
            scatter(x, auc, 50, 'b*');
            xlim([0, N+1]); ylim([0.5, 1]);
            ylabel('Area Under Curve')
            text(x-0.3,auc,num2str(auc, '%.2f'), 'FontSize', 14);
        end

        yyaxis right
        scatter(x, accur, 50 ,'ro');
        ylabel('Accuracy %')
        text(x+0.1,accur,num2str(accur, '%.2f'), 'FontSize', 14);
        ylim([50, 100]);

        hold off
        title(figTitle)
        xlabel('Input Dataset')
        xticks(x); xticklabels(lineNames); xtickangle(45); 
        set(gcf, 'Position', get(0, 'Screensize'));

    case 'visual'
%% Visualization of malignancy score %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cmapSize = 100; % default size of 60 shows visible discretization
        if ischar(cmap)

            try
                cmap = eval([cmap '(' num2str(cmapSize) ');']);
            catch
                fprintf('Colormap ''%s'' is not supported. Using ''jet''.\n',cmapName);
                cmap = jet(cmapSize);
            end
        end
        colormap(cmap);

        clf(gcf);
        axes(gcf);

        hB = imagesc(I);axis image off;
        climF = [min(overlay(:)), max(overlay(:))];
        
        if ~strcmp(figTitle, 'Ground Truth')
            % Add the front image on top of the back image
            hold on;
            hF = imagesc(overlay); %, climF);

            % If images are different sizes, map the front image to back coordinates
            set(hF,'XData',get(hB,'XData'),...
                'YData',get(hB,'YData'))

            % Make the foreground image transparent
            alphadata = alpha.*(overlay > climF(1));
            set(hF,'AlphaData',alphadata);

            c = colorbar('location','southoutside', 'Ticks', [0 0.5 1], 'TickLabels', {'low', 'medium', 'high'});
            c.Label.String = 'Malignancy Probability';
            c.Label.FontSize = 10;
            c.Label.FontWeight = 'bold';
            c.Limits = [0,1];
            set(gcf,'Visible','on');
        end
        title(figTitle)
        if exist('coordinates', 'var') && ~isempty(coordinates)
            hold on
            for i = 1:size(coordinates, 1)
                x = coordinates(i,1);
                y = coordinates(i,2);
                if (curves(i) == 1) 
                    plot(x, y,'rx', 'MarkerSize', 10, 'MarkerEdgeColor', 'r', 'LineWidth', 3);
                else 
                    plot(x, y,'go', 'MarkerSize', 10, 'MarkerEdgeColor', 'g', 'LineWidth', 3);
                end
            end
            hold off
        end
        
    case 'classificationPerformance'
%% Plot classification performance comparison %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        grouping = lineNames;
        aucVal = performance(1,:)';
        accVal = performance(2,:)';
        if (size(grouping,2) > 1) 
            gscatter(aucVal, accVal, grouping, 'rrrrcccc','.*^x+', [10, 10]);
        else 
            gscatter(aucVal, accVal, grouping, 'rbgck', '.x', 10);
        end
        xlabel('AUC');
        ylabel('Accuracy');
        title('Perfomance of Top Classifiers', 'FontSize', 15)
        legend('Location','northeastoutside')
        ylim([0.7, 0.9])
        xlim([0.7, 0.9])

    case 'featureImportance'
%% Plot important features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         features = categorical({'760', '635', '620', '650', '645', '655', '665', '625', '660'});
         importance = [0.02170081 0.03423029 0.04734384 0.07019354 0.08864051 0.1096366 0.13921365 0.14391144 0.17794282];
         bar(features, importance);
         title('Top 10 Important Wavelengths', 'FontSize', 15);
         ylabel('Relative Importance', 'FontSize', 12);
         xlabel('Wavelength (nm)', 'FontSize', 12);
         plotName = fullfile(saveOptions.savedir,'Importance', 'featimp');

    case {'classificationPerformanceBars' }
%% Plot bar performance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if (saveInBW)
            color(1,:)=HInt2RGB(1,100); % red, darkest 
            color(2,:)=HInt2RGB(3,64); % green, less dark 
            color(3,:)=HInt2RGB(7,10); % blue cyan, lightest
        else 
            color = hsv(3);
            %color = flip(color);
        end
        
        plotN = length(figTitle)-1;
        superTitle = figTitle{1};
        figTitle = figTitle(2:end); 
        for ploti=1:plotN 
            subplot(plotN,1,ploti);
            figTitleCur = figTitle{ploti};
            performanceCur = performance{ploti};
            hasLegend = ploti == 1;
            if contains(figTitleCur, 'classifier')
                legTitle = 'Classifier';
                figTitleCur = 'Cross Validated AUC';
                labelx = 'Feature Vector';
                labely = 'ROC AUC';
                cs = {'SVM', 'KNN', 'RF'};
                lims = [70, 100];
            elseif contains(figTitleCur, 'fixing')
                legTitle = 'Tissue';
                figTitleCur = 'Cross Validated AUC';
                labelx = 'Feature Vector';
                labely = 'ROC AUC';
                cs = {'Unfixed', 'Fixed', 'Mixed'};
                lims = [70, 100];
            else 
                legTitle = 'Tissue';
                cs = {'Unfixed', 'Fixed', 'Mixed'};
                labelx = 'Classifier';
                labely = strcat([figTitleCur, '(%)']);
                figTitleCur = '';
                %figTitleCur = strcat(['Cross Validated ', figTitleCur]);
                lims = [50, 100];
            end
            barPlot = GetBarPlot(lineNames, performanceCur * 100, color, figTitleCur, labelx, labely, cs, legTitle, lims, 9, hasLegend);
        end 
        suptitle(superTitle);
        plotName = fullfile(saveOptions.savedir, superTitle);

        
    case 'reconstructionPerformanceBars'
%% Plot bar performance for reconstruction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if (saveInBW)
            color(1,:)=HInt2RGB(1,100); % red, darkest 
            color(2,:)=HInt2RGB(3,64); % green, less dark 
            color(3,:)=HInt2RGB(7,10); % blue cyan, lightest
        else 
            color = hsv(3);
            %color = flip(color);
        end
        plotName = fullfile(saveOptions.savedir, 'bar_plots', 'bar_reconstruction');
        barPlot = GetBarPlot(lineNames, performance * 100, color, figTitle, 'Malignancy', 'NRMSE(%)',...
            {'unfixed', 'fixed', 'sectioned'}, {'Fixing'}, [0, 15], 12, true);
            
    otherwise 
        error('Unsupported plot');
end

%% SaveImages%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get rid of tool bar and pulldown menus that are along top of figure.
%set(gcf, 'Toolbar', 'none', 'Menu', 'none');

%set(0, 'CurrentFigure', fig);
%handle = get(groot,'CurrentFigure');

if (savePlot && ~isempty(plotName))
    filename = strrep(plotName, '.mat', '');
    
    [filepath,name,~] = fileparts(filename);
    filepathBW = fullfile(filepath, 'bw');    
    mkdir_custom(filepath);
    mkdir_custom(filepathBW);

    
    if (saveInHQ)
        filename = fullfile(filepath, strcat(name, '.png'));
        export_fig(filename , '-png','-native', '-nocrop');
        %print(handle, strcat(plotName, '.png'), '-dpng', '-r600');
    else
        filename = fullfile(filepath, strcat(name, '.png'));
        saveas(fig , filename);
        filename = fullfile(filepath, strcat(name, '.eps'));
        export_fig(filename , '-eps', '-transparent', '-r900',  '-RGB');
        if ~(saveInBW)
            filename = fullfile(filepathBW, strcat(name, '.eps'));
            export_fig(filename , '-eps', '-transparent', '-r900',  '-gray');
        else
            plots(plotType, fig, [], '', 'Wavelength', wavelength, 'Illumination', illumination, 'Sensitivity', sensitivity, 'SaveOptions', saveOptions, 'saveInBW', true);
        end
    end
end
warning(w);

end

function [outname] = strrepAll(inname)
    
    [~, outname] = fileparts(inname);
    outname = strrep(outname, '\', ' ');
    outname = strrep(outname, '_', ' ');
    outname = strrep(outname, '.csv', '');
    outname = strrep(outname, '.mat', '');

end

function [sampleName] = simpleSampleName(name)
    attr = strsplit(name);
    sampleName = attr{1};
end

function [barPlot] = GetBarPlot(lineNames, values, color, figTitle, labelx, labely, categories, legendTitle, ylims, ticklabelsize, hasLegend)

            barPlot = bar(categorical(lineNames), values, 'FaceColor', 'flat');
            for k = 1:length(categories)
                barPlot(k).CData = color(k,:);
            end
            hleg = legend(categories, 'FontSize', 9.5, 'Location', 'eastoutside');
            htitle = get(hleg,'Title');
            set(htitle,'String',legendTitle)
            if ~(hasLegend)
                hleg.Visible = 'off';
            end
            ax = gca;
            ax.FontSize = ticklabelsize; 
            ylim(ylims);
            xlabel(labelx, 'FontSize', 12);
            ylabel(labely, 'FontSize', 12);
            if ~strcmp(figTitle, '')
                title(figTitle, 'FontSize', 16); % should be after gca, set to work
            end
end
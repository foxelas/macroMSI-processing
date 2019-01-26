function [] = plots(plotType, varargin)
%PLOTS Each time plots an informative plot based on 'plotType'
%
%Syntax:
%   plots(plotType, 1, curves, 'leftcamera23')
%
%Available plots:
%'sensitivity',
%'illuminationAndSensitivity'
%'illumination'
%'estimationComparison'
%'allEstimations'
%'singlemeasurement'
%'d65'
%'methodErrors'

%%%%%%%%%%%%%%%%%%%%%%%%Initialize%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
defaultfc = [450, 465, 505, 525, 575, 605, 630];
addRequired(p, 'plotType', @(x) any(validatestring(x, {'sensitivity', 'illuminationAndSensitivity', 'illumination', 'estimationComparison', ...
    'allEstimations', 'singlemeasurement', 'd65', 'methodErrors', 'overlapSpectrum', 'overlapSpectrumSample', 'pca', 'lda', 'pca b', 'lda b', ...
    'pcalda', 'classificationErrors', 'cropped', 'segmentation'})));
addOptional(p, 'fig', -1);
addOptional(p, 'curves', []);
addOptional(p, 'name', '', @(x) ischar(x));
addParameter(p, 'Wavelength', []);
addParameter(p, 'Sensitivity', []);
addParameter(p, 'Illumination', []);
addParameter(p, 'Method', '', @(x) ischar(x));
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

parse(p, plotType, varargin{:});
plotType = p.Results.plotType;
fig = p.Results.fig;
curves = p.Results.curves;
name = p.Results.name;
sensitivity = p.Results.Sensitivity;
wavelength = p.Results.Wavelength;
illumination = p.Results.Illumination;
method = p.Results.Method;
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
end

if ~isempty(lineNames)
    lineNames = cellfun(@(x) strrep(x, '_', ' '), lineNames, 'un', 0);
end
if ~isempty(name)
    name = strrep(name, '_', ' ');
end
if ~isempty(method)
    method = strrep(method, '_', ' ');
end

w = warning('off', 'all');
switch plotType
    case 'sensitivity'     
        %% Plot camera sensitivity%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sensitivityN = size(sensitivity, 2);
        hold on;
        for j = 1:sensitivityN
            plot(wavelength, sensitivity(:, j));
        end
        hold off;
        title('Camera Sensitivities');
        xlabel('Wavelength \lambda (nm)');
        ylabel('Sensitivity');
        xlim([400, 700]);
        if (sensitivityN == 7)
            fnames = cellstr(num2str(defaultfc'));
        else
            fnames = {'red channel'; 'green channel'; 'blue channel'};
        end
        legend(fnames)
        % End Plot camera sensitivity
        
    case 'illuminationAndSensitivity'     
        %% Plot illumination and sensitivity%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %{
            % Arguments  ( wavelength, bandWavelength, illumination )
            % wavelength: 401x1 or 81x1
            % bandWavelength: the names of the bands
            % illumination: the illumination light 401x7
        %}
        cl = colormap(lines);
        yyaxis left;
        hold on;
        for j = 1:numel(fc)
            plot(wavelength, illumination(:, j), 'DisplayName', ['illum', num2str(fc(j)), ' nm'], 'LineStyle', '-', 'Marker', 'none', 'Color', cl(j, :));
        end
        hold off;
        title('Luminous intensity');
        xlabel('Wavelength \lambda (nm)');
        ylabel('Luminous Intensity (cd/m^2)');
        
        yyaxis right;
        sensitivityN = size(sensitivity, 2);
        if (sensitivityN == 7)
            fnames = cellstr(num2str(defaultfc'));
        else
            fnames = {'red channel'; 'green channel'; 'blue channel'};
        end
        hold on;
        for j = 1:sensitivityN
            plot(wavelength, sensitivity(:, j), 'DisplayName', ['sens', fnames{j}, 'nm'], 'LineWidth', 2, 'LineStyle', '-', 'Marker', 'none', 'Color', cl(j+numel(fc), :));
        end
        hold off;
        ylabel('Camera sensitivity');
        xlim([400, 700]);
        title('Overlap of illumination and camera sensitivity spectrum')
        legend;
        
        %End Plot illumination and sensitivity
        
    case 'illumination'
        %% Plot illumination%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                     
        % Arguments ( wavelength, fc, illumination )
        % wavelength: 401x1 or 81x1
        % fc: central frequencies of the bands (string cell )
        % illumination: the illumination light 401x7
        hold on;
        for j = 1:numel(fc)
            plot(wavelength, illumination(:, j), 'DisplayName', [num2str(fc(j)), ' nm']);
        end
        hold off;
        legend;
        xlim([400, 700]);
        title('Luminous intensity');
        xlabel('Wavelength \lambda (nm)');
        ylabel('Luminous Intensity (cd/m^2)');
        
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
        end
        
        if isempty(markers)
            for i = 1:curveN
                markers{i} = 'none';
            end
        end
        
        if ~exist('MSIreflectances', 'var')
            plotN = 2;
        else
            plotN = 3;
        end
        
        % Set up figure properties:
        % Enlarge figure to full screen.
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 1, 1]);
        % Get rid of tool bar and pulldown menus that are along top of figure.
        set(gcf, 'Toolbar', 'none', 'Menu', 'none');
        
        figTitle = sprintf('Sample: %s | Method: %s',name, method);
        %marker = {'none', 'o', '+', '*', '.', 'none', 'o', '+', '*', '.'};
        color = colorcube(curveN+10); %color = [ 'b', 'g' , 'k',  'c' , 'y' , 'm'];
        
        [~, peakIdx] = ismember(fc, wavelength); % mark filter wavelengths.
        subplot(1, plotN, [1, 2]);
        hold on
        for i = 1:curveN
            plot(wavelength, curves(:, i), 'Color', color(i, :), 'Marker', markers{i}, 'LineWidth', 1.3, 'DisplayName', lineNames{i + 1}); % plot estimated reflectances
        end
        plot(wavelength(peakIdx), curves(peakIdx, 1), 'rx', 'DisplayName', lineNames{1}, 'LineWidth', 1.3); % plot measured reflectance
        hold off
        
        xlabel('Wavelength \lambda (nm)');
        ylabel('Reflectance Spectrum');
        xlim([400, 700]);
        suptitle('Comparative plot of Wiener estimation results')
        title(figTitle);
        legend({lineNames{2:(curveN+1)}, lineNames{1}}, 'Location', 'northwest', 'FontSize', 14); % 'Orientation','horizontal');
        
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
            sprintf('Sample %s with %s smoothing matrix ', name, method)};
        figTitle = strrepAll(figTitle);
        title(figTitle);
        xlabel('Wavelength \lambda (nm)');
        ylabel('Reflectance');
        legend(h1, {'Measured reflectance'});
        legend(h2, {'Estimated reflectance'});
        % end plots all estimated curves for all pixels in an image area
        
    case 'singlemeasurement'       
        %% plot single measurement %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        hold on
        plot(wavelength, curves, 'DisplayName', strrep(lineNames, '_', ' '), 'Color', 'm')
        xlabel('wavelength (nm)');
        ylabel('Reflectance ratio');
        title('Reflectance spectrum of a point object');
        xlim([400, 700]);
        hold off
        %end plot single measurement
        
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
        if contains(name, 'nmse', 'IgnoreCase', true)
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
        
        if contains(name, 'System') || contains(name, 'Matrix')
            smms = unique({errors.options.smoothingMatrixMethod}, 'stable');
            pvsms = unique({errors.options.pixelValueSelectionMethod}, 'stable');
            if (length(smms) == 1)
                [pvsms, smms] = deal(smms, pvsms);
                labelx = 'Pixel value selection method';
            else
                labelx = 'Wiener estimation variation';
            end
            
        elseif contains(name, 'Noise')
            smms = unique({errors.options.noiseType}, 'stable');
            pvsms = unique({errors.options.pixelValueSelectionMethod}, 'stable');
            labelx = 'Noise model';
            
        end
        
        marker = ['o', 's', 'd', '^', '*', 'h', 'p', 'v', '<', '+', '>'];
        if ((length(smms) > 1) && (length(pvsms) > 1))
            pvsmColor = colorcube(length(pvsms)+10);
        else
            pvsmColor = jet(length(pvsms)+10);
        end
        
        hold on
        % Generate dummy info for plot handles "h"
        h = zeros(length(pvsms), 1);
        for i = 1:length(pvsms)
            h(i) = plot([NaN, NaN], 'Color', pvsmColor(i, :), 'Marker', marker(i), 'DisplayName', pvsms{i});
        end
        if contains(name, 'avg', 'IgnoreCase', true)
            i = i + 1;
            h(i) = plot([NaN, NaN], 'ko', 'DisplayName', 'Mean and Std');
        end
        if contains(name, 'max', 'IgnoreCase', true)
            i = i + 1;
            h(i) = plot([NaN, NaN], 'k-.', 'DisplayName', 'Max value');
        end
        if contains(name, 'min', 'IgnoreCase', true)
            i = i + 1;
            h(i) = plot([NaN, NaN], 'k:', 'DisplayName', 'Min value');
        end
        % end dummy data for legend entries
        
        for i = 1:length(pvsms)
            for j = 1:length(smms)
                if contains(name, 'avg', 'IgnoreCase', true)
                    errorbar(j, avge((i - 1)*length(smms)+j), stde((i - 1)*length(smms)+j), 'Color', pvsmColor(i, :), 'LineWidth', 1, 'Marker', marker(i), 'MarkerSize', 10);
                end
            end
            if contains(name, 'max', 'IgnoreCase', true)
                plot(1:length(smms), maxe, 'Color', pvsmColor(i, :), 'LineStyle', '-.');
            end
            if contains(name, 'min', 'IgnoreCase', true)
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
        legend(h, 'Location', 'best');
        title('Comparison of estimation results for the various method variations')
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
        
        
    case {'pca', 'lda', 'pca b', 'lda b', 'pcalda'}     
        %% plot discriminant analysis results%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if strcmp(plotType, 'pca')
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
        attr = split(lineNames, '_');
        sample = {attr{:, 1}}';
        type = {attr{:, 2}}';
        isNormal = {attr{:, 3}}';
        idx = strcmp(isNormal, 'Normal') | strcmp(isNormal, 'Benign');
        colors(idx) = 'b';
        colors(~ismember(idx, 1:length(isNormal))) = 'r';
        
        %dummy legends
        h = zeros(2, 1);
        hold on
        h(1) = plot([NaN, NaN], 'Color', 'r', 'DisplayName', 'Cancer');
        h(2) = plot([NaN, NaN], 'Color', 'b', 'DisplayName', 'Normal');
        hold off
        
        if contains(name, 'Fix', 'IgnoreCase', true)
            hold on
            h(3) = plot([NaN, NaN], 'Color', 'k', 'Marker', marker(1), 'DisplayName', 'Fixed');
            h(4) = plot([NaN, NaN], 'Color', 'k', 'Marker', marker(2), 'DisplayName', 'Unfixed');
            hold off
            markers = repmat(marker(1), observations, 1);
            idx = arrayfun(@(x) any(strcmp(x, 'unfixed')), type);
            markers(idx) = marker(2);
        end
        
        if contains(name, 'Sample', 'IgnoreCase', true)
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
        titl = strsplit(plotName, '\');
        figTitle = strrepAll(titl{end});
        suptitle(figTitle);
        
        legend(h, 'Location', 'best');
        set(gcf, 'Position', get(0, 'Screensize'));
        % End plot discriminant analysis results 
        
    case 'cropped'
        %% Show cropped section %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        markers = {'r*', 'g*', 'm*', 'y*', 'c*'};
        imshow(I);
        hold on
        for i = 1:size(coordinates, 1)
            x = coordinates(i,1);
            y = coordinates(i,2);
            plot(x, y, markers{i}, 'LineWidth', 2, 'MarkerSize', 5);
        end
        hold off
        figTitle = strrepAll(strcat('Cropped area of ', plotName));
        title(figTitle);
        
    case 'segmentation'
        %% Show segmentation images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Colors defined from https://academo.org/demos/wavelength-to-colour-relationship/
        bandColors = [  0,70,255   ; 
                0,146,255  ;
                0,255,84   ;
                74,255,0   ;
                240,255,0  ;
                255,173,0  ;
                255,79,0   ;
                255,255,255
              ];
        bandColors = bandColors./255;
        colormap(bandColors);
        imshow(I);
        c = colorbar('location','southoutside', 'Ticks', linspace(0,1,9),...
             'TickLabels',{'450','465','505','525','575','605', '630', 'All', ''});
        c.Label.String = 'Respective MSI band (nm)';
        figTitle = strrepAll(strcat('Segmented area of ', plotName));
        title(figTitle);
        plotName = strcat(plotName, '_segments');
end
warning(w);

%% SaveImages%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (savePlot && ~isempty(plotName))
    plotName = strrep(plotName, '.mat', '');
    % plotName = strrep(plotName, ' ', ''); %         plotName = strrep(plotName, '_', '');
    strIdxs = strfind(plotName, '\');
    strIdx = strIdxs(end);
    fn = plotName(1:strIdx);
    
    if ~exist(fn, 'dir')
        mkdir(fn);
        addpath(fn);
    end
    
    set(0, 'CurrentFigure', fig);
    if (saveInHQ)
        export_fig(strcat(plotName, '.jpg') , '-jpg','-native');
        %print(fig, strcat(plotName, '.jpg'), '-djpeg', '-r600');
    else
        export_fig(strcat(plotName, '.jpg') , '-jpg');
        %print(fig, strcat(plotName, '.jpg'), '-djpeg');
    end
end

end

function [outname] = strrepAll(inname)
    
    [~, outname] = fileparts(inname);
    outname = strrep(outname, '\', ' ');
    outname = strrep(outname, '_', ' ');
    outname = strrep(outname, '.csv', '');
    outname = strrep(outname, '.mat', '');

end

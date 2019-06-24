%Start of reflectance estimation comparison
options.action = 'Refest_Preset_plusrgb';

%% Comparison settings
[smms, pvsms, nms, plotType] = getComparisonSets(options);
plusRGBAnalysis =  contains(lower(options.action), 'plusrgb');
hasRGBAnalysis = strcmp(options.pixelValueSelectionMethod, 'rgb') || contains(lower(options.action), 'rgb');
isNoComparison = (contains(lower(options.action), 'preset') || contains(lower(options.action), 'simple'));
[lineNames, methodsN, pixelValueSelectionMethods, smoothingMatrixMethods, noiseTypes, smmsN, pvsmsN, nmsN] = getLineNames(pvsms, smms, nms);

%% Comparison
w = warning('off', 'all');

tic;
estimatedSpectra = zeros(size(Spectra));
estimatedRGBSpectra = zeros(size(Spectra));
errFixingInfo = struct('RMSE', [], 'NMSE', [], 'fixing', [], 'malignancy', []);
errFixingInfoRGB = struct('RMSE', [], 'NMSE', [], 'fixing', [], 'malignancy', []);

for k = 1:msiN
    % Retrieve MSI data
    
    infile = fullfile(options.systemdir, 'infiles', strcat('poi_', num2str(k), '.mat'));
    load(infile, 'poiName', 'poiRAW', 'poiSegmentMask', ...
        'roiSeeds', 'measuredSpectrum', 'poiWhite');
        
    msi = poiRAW;
    mask = poiSegmentMask;
    measured = measuredSpectrum;
    if (hasRGBAnalysis); rgb = readRGBImage(poiWhite); end
        
    % Estimation
    estimated = zeros(methodsN, length(wavelength));
    
    for j = 1:methodsN
        [options, optionsRgb]= updateOptions(options, j, pvsms, smms, nms);
        
        if (hasRGBAnalysis) %Case with only RGB
            [estRGB, rmseRGB, nmseRGB, ~] = reflectanceEstimation(rgb, mask, measured, ID(k), optionsRgb);
            estimatedRGBSpectra(k,:) = estRGB;
            errFixingInfoRGB(k) = struct('RMSE', rmseRGB, 'NMSE', nmseRGB, 'fixing', ID(k).Type, 'malignancy', ~ID(k).IsBenign);

        end
        [estimated(j,:), rmse, nmse, ~] = reflectanceEstimation(msi, mask, measured, ID(k), options);                 
        
        if (isNoComparison)
            estimatedSpectra(k,:) = estimated(j,:);
            errFixingInfo(k) = struct('RMSE', rmse, 'NMSE', nmse, 'fixing', ID(k).Type, 'malignancy', ~ID(k).IsBenign);
        end
      
    end
    
    if (options.showImages)       
        if (~plusRGBAnalysis) ; lines = [measured; estimated]; else;  lines = [measured; estimated; estRGB]; end   
        options.saveOptions.plotName = fullfile(options.saveOptions.savedir, '5-ReflectanceEstimation', options.action,...
        strcat( options.action, '_', num2str(ID(k).Index)));  
        plotReconstructedCurves(lines', lineNames, wavelength, 'Reflectance Estimation Comparison',...
            1,options.saveOptions);        
        pause(0.1)
    end
    
end
warning(w);

%% Export results
errors = GetErrorInfoStruct(rmse, nmse, pixelValueSelectionMethods, smoothingMatrixMethods, noiseTypes);

if (isNoComparison)
    fprintf('%s\n', options.noiseType);
    % show NRMSE error per set
    options.action = strjoin({options.action, options.noiseType}, '_');
    if strcmp(options.smoothingMatrixMethod, 'adaptive')
        options.action = strjoin({options.action, 'adaptive'}, '_');
    end
    nrmse = showNRMSE([errFixingInfo.RMSE], Spectra, ID);
    titl = 'MSI-Reconstruction NRMSE';
    options.saveOptions.plotName = fullfile(options.saveOptions.savedir, '6-ReflectanceEstimationPerformance', options.action, titl);
    plotReconstructionPerformanceBars(nrmse,{'benign', 'malignant'},titl,2,options.saveOptions);
    writeOutput(options, errFixingInfo, estimatedSpectra, false);
    nrmse = showNRMSE([errFixingInfoRGB.RMSE], Spectra, ID);
    
    titl = 'RGB-Reconstruction NRMSE';
    options.saveOptions.plotName = fullfile(options.saveOptions.savedir, '6-ReflectanceEstimationPerformance', options.action, titl);
    plotReconstructionPerformanceBars(nrmse,{'benign', 'malignant'},titl,3,options.saveOptions);
    writeOutput(options, errFixingInfoRGB, estimatedRGBSpectra, true);

else
    writetable(struct2table(errors), strcat('../../logs/reflectance_estimation_comparison_', options.dataset,'_log.xlsx'));

    options.saveOptions.plotName = generateName(['ComparisonRmse', options, plotType]);
    plots('methodErrors', 3, [], ['Rmse', plotType], 'Errors', errors, 'SaveOptions', options.saveOptions)

    options.saveOptions.plotName = generateName(['ComparisonNmse', plotType], options);
    plots('methodErrors', 4, [], ['Nmse', plotType], 'Errors', errors, 'SaveOptions', options.saveOptions)
    
    save(generateName(['ComparisonErrors', plotType, '.mat'], options), 'errors');
end
t2 = toc; 
fprintf('Run time time %.5f\n', t2);

%End of reflectance estimation comparison


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function errorStrutct = GetErrorInfoStruct(rmse, nmse, pixelValueSelectionMethods, smoothingMatrixMethods, noiseTypes)
    errorStrutct = struct('avgrmse', num2cell(mean(rmse, 2)), 'minrmse', num2cell(min(rmse, [], 2)), 'maxrmse', num2cell(max(rmse, [], 2)), ...
                'stdrmse', num2cell(std(rmse, [], 2)), 'avgnmse', num2cell(mean(nmse, 2)), 'minnmse', num2cell(min(nmse, [], 2)), ...
                'maxnmse', num2cell(max(nmse, [], 2)), 'stdnmse', num2cell(std(nmse, [], 2)), 'stdeSEM',  num2cell(std(rmse, [], 2)./size(rmse,2)), ...
                'pixelValueSelectionMethod', pixelValueSelectionMethods, ...
                'smoothingMatrixMethod', smoothingMatrixMethods, 'noiseType', noiseTypes);
end

function [lineNames, methodsN, pixelValueSelectionMethods, smoothingMatrixMethods, noiseTypes, smmsN, pvsmsN, nmsN]= getLineNames(pvsms, smms, nms)
    smmsN = numel(smms);
    pvsmsN = numel(pvsms);
    nmsN = length(nms);
    methodsN = smmsN * pvsmsN * nmsN; 
    
    pixelValueSelectionMethods = cell(methodsN,1);
    smoothingMatrixMethods = cell(methodsN,1);
    noiseTypes = cell(methodsN,1);

    lineNames = {'Channel \lambda', 'Measured'};
    for l = 1:pvsmsN
        for n = 1:smmsN
            for m = 1:nmsN
                j = sub2ind([ pvsmsN , smmsN, nmsN],  l, n, m);  
                pixelValueSelectionMethods{j} = pvsms{l};
                smoothingMatrixMethods{j} = smms{n};
                noiseTypes{j} = nms{m};  
                if (pvsmsN == 1) && (smmsN == 1) && (nmsN == 1) 
                    lineNames{j+2} = strcat('MSI Estimate');
                elseif (pvsmsN < 2)
                    lineNames{j+2} = strcat('MSI Estimate -', smms{n});
                elseif (smmsN < 2)
                    lineNames{j+2} = strcat('MSI Estimate -', pvsms{l});
                else 
                    lineNames{j+2} = strcat('MSI Estimate -', pvsms{l}, '|', smms{n});
                end
                if (nmsN > 1)
                    lineNames{j+2} = strcat('MSI Estimate -', nms{m});
                end
            end
        end
    end
    lineNames{end + 1} = 'RGB estimate';
end

function [updatedOptions, optionsRgb] = updateOptions(options, j, pvsms, smms, nms)
    smmsN = numel(smms);
    pvsmsN = numel(pvsms);
    nmsN = length(nms);
    [l,n,m] = ind2sub([ pvsmsN , smmsN, nmsN], j);
    options.pixelValueSelectionMethod = pvsms{l};
    options.smoothingMatrixMethod = smms{n};
    options.noiseType = nms{m};   
    optionsRgb = options;
    updatedOptions = options;

    % In the case of 3channel RGB image
    if strcmp(options.pixelValueSelectionMethod, 'rgb') || contains(lower(options.action), 'rgb')
        optionsRgb.pixelValueSelectionMethod = 'rgb';
        optionsRgb.noiseType = 'fromOlympus';
    end
end

function [rgb] = readRGBImage(whiteImg)
    load('saved parameters\color_correction.mat', 'illuminant_gw1');
    tempRGB = chromadapt(whiteImg, illuminant_gw1, 'ColorSpace', 'linear-rgb'); %color adjustment
    rgb = reshape(tempRGB, [3, size(tempRGB, 1), size(tempRGB, 2)]);
end

function [smms, pvsms, nms, plotType] = getComparisonSets(options)

    smmsFull = { 'Cor_All', 'Cor_Fixing', 'Cor_Sample', 'adaptive'}; %, 'Cor_SampleMalignancy', 'Cor_SampleMalignancyFixing', 'Cor_MalignancyFixing', 'Cor_SampleMalignancyFixing', 'markovian', 'Cor_Macbeth',};
    pvsmsFull = {'green', 'rms', 'adjusted', 'extended', 'rgb'};
    nmsFull = {'sameForChannel 0.0001', 'sameForChannel 0.001', 'sameForChannel 0.01',  'sameForChannel 0.1', 'sameForChannel 0.5',...
        'diffForChannel 0.0016, 0.0016, 0.0008, 0.0005, 0.0008, 0.0148, 0.0015', 'diffForChannel 0.16, 0.16, 0.08, 0.05, 0.08, 1.48, 0.15',...
        'SNR 25', 'SNR 20', 'SNR 15', 'SNR 5',  ...
        'spatial 0.002 0.04', 'spatial 0.005 0.5', 'spatial 0.002 0.01', 'fromOlympus'}; %, 'fromOlympus', 'spatial', 'spatial 0.015 0.015' 

    if contains(lower(options.action), 'matrixsystem')
            smms = smmsFull;
            pvsms = pvsmsFull;
            plotType = 'MatrixSystemAvg';        
            nms = {'SNR'}; 

    elseif contains(lower(options.action), 'matrixnoise')
            smms = smmsFull;
            pvsms = {'extended'};
            plotType = 'MatrixNoiseAvg';        
            nms = nmsFull; 

    elseif contains(lower(options.action), 'smoothing matrix') || contains(lower(options.action), 'matrix')
            smms = smmsFull;
            pvsms = {'extended'};
            plotType = 'MatrixAvgMinMax';
            nms = {'SNR'};

    elseif contains(lower(options.action), 'color system') || contains(lower(options.action), 'system')   
            smms = {'Cor_Malignancy'};
            pvsms = pvsmsFull;
            plotType = 'SystemAvgMinMax';
            nms = {'SNR'};

    elseif contains(lower(options.action), 'noise model') || contains(lower(options.action), 'noise')   
            pvsms = {'extended'};
            smms = {'Cor_SampleMalignancy'};
            plotType = 'NoiseAvgMinMax'; %Or 'NoiseAvgMinMaxWithNone'
            nms = nmsFull;

    elseif contains(lower(options.action), 'simple') 
            pvsms = {'extended'};
            smms = {options.smoothingMatrixMethod};
            nms = {options.noiseType};
            plotType = '';

    elseif contains(lower(options.action), 'preset')
            pvsms = {'extended'}; %{'extended'};
            smms = {'Cor_Sample'}; %{'adaptive'}; % {'Cor_Sample'}; %{'Cor_Malignancy'};  {'Cor_All'}
            nms =  {'sameForChannel 0.1'}; %{'sameForChannel 0.1'}; %{'fromOlympus'};%{'spatial olympus'};
            plotType = '';

    elseif contains(lower(options.action), 'extra') 
            pvsms = {'extended'};
            smms = {'Cor_All', 'adaptive'};
            nms = {'sameForChannel 0.0001', 'spatial 0.0004 0.002'};
            plotType = '';

    else 
            error('Not implemented yet.')
    end

end

function nrmseArray = showNRMSE(rmse, measured, ID)
    i = 0;
    nrmseArray = zeros(6, 1);
    for state = {'unfixed', 'fixed', 'cut'}
        for malignancy = 0:1
            i = i + 1;
            nrmseArray(i) = GetNRMSE(measured, rmse, [ID.IsBenign] ~= malignancy & strcmp({ID.Type}, state));
            %fprintf('NRMSE for malignancy %d  and state %s = %.4f\n',malignancy, state{1}, nrmseArray(i));               
        end
        nrmse = GetNRMSE(measured, rmse, strcmp({ID.Type}, state));
        %fprintf('NRMSE for %s = %.4f\n', state{1}, nrmse);
    end
   nrmseArray = reshape(nrmseArray,[2, 3]); 
   
    for malignancy = 0:1
        nrmse = GetNRMSE(measured, rmse, [ID.IsBenign] ~= malignancy);
        %fprintf('NRMSE overall for malignancy %d = %.4f\n',malignancy, nrmse);
    end
    
    nrmse = mean(rmse) /  ( max(measured(:)) - min(measured(:)));
    fprintf('NRMSE overall = %.4f\n', nrmse);
end

function nrmse = GetNRMSE(measured, rmse, stateId)
    currState = measured(stateId, :);
    maxDif = max(currState(:)) - min(currState(:));
    nrmse = mean(rmse(stateId)) / maxDif;
end

function [] = writeOutput(options, errFixingInfo, estimatedSpectra, hasRGBAnalysis)
    
    if (hasRGBAnalysis)
        writetable(struct2table(errFixingInfo), strcat('../../logs/reflectance_estimation_comparison_', options.dataset,'_log.xlsx'), 'Sheet', 'RGB');
        EstimatedRGBSpectra = estimatedSpectra;
    else
        writetable(struct2table(errFixingInfo), strcat('../../logs/reflectance_estimation_comparison_', options.dataset,'_log.xlsx'), 'Sheet', 'MSI');
        EstimatedSpectra = estimatedSpectra;
    end
    filename = mkdir_custom(fullfile(options.saveOptions.savedir, '5-ReflectanceEstimation', options.action, 'refest.mat'));
    if exist(filename, 'file')
        if (hasRGBAnalysis)
            save(filename, 'EstimatedRGBSpectra','-append');
        else
            save(filename, 'EstimatedSpectra', '-append');
        end
    else
        if (hasRGBAnalysis)
            save(filename, 'EstimatedRGBSpectra');
        else
            save(filename, 'EstimatedSpectra');
        end
    end
end

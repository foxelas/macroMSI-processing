%Start of reflectance estimation comparison
%% Comparison settings
[smms, pvsms, nms, plotType] = getComparisonSets(options);
plusRGBAnalysis =  contains(lower(options.action), 'plusrgb');
hasRGBAnalysis = strcmp(options.pixelValueSelectionMethod, 'rgb') || contains(lower(options.action), 'rgb');
isNoComparison = (contains(lower(options.action), 'preset') || contains(lower(options.action), 'simple')) && ~contains(lower(options.action), 'plusrgb');
[lineNames, methodsN, pixelValueSelectionMethods, smoothingMatrixMethods, noiseTypes, smmsN, pvsmsN, nmsN] = getLineNames(pvsms, smms, nms);

%% Comparison
w = warning('off', 'all');

tic;
rmse = zeros(methodsN, msiN);
nmse = zeros(methodsN, msiN);
estimatedSpectra = zeros(size(Spectra));
errFixingInfo = struct('rmse', [], 'fixing', [], 'malignancy', []);
for k = 1:msiN
    % Retrieve MSI data
    msi = MSIs{k};
    mask = Masks{k};
    measured = Spectra(k,:);
    if (hasRGBAnalysis); rgb = readRGBImage(WhiteIs{k}); end
        
    % Estimation
    estimated = zeros(methodsN, length(wavelength));
    
    for j = 1:methodsN
        [options, optionsRgb]= updateOptions(options, j, pvsms, smms, nms);
        
        if (hasRGBAnalysis && ~plusRGBAnalysis) %Case with only RGB
            [estimated(j,:), rmse(j, k), nmse(j, k), ~] = reflectanceEstimation(rgb, mask, measured, ID(k), optionsRgb);

        elseif (~hasRGBAnalysis) %Case with only MSI 
            [estimated(j,:), rmse(j, k), nmse(j, k), ~] = reflectanceEstimation(msi, mask, measured, ID(k), options);                 

        else % Case with both RGB and MSI 
            [estRGB, ~, ~, ~] = reflectanceEstimation(rgb, mask, measured, ID(k), optionsRgb);
            [estimated(j,:), rmse(j, k), nmse(j, k), ~] = reflectanceEstimation(msi, mask, measured, ID(k), options);                 
        end
        
        if (isNoComparison)
            estimatedSpectra(k,:) = estimated(j,:);
            errFixingInfo(k) = struct('rmse', rmse(j, k), 'fixing', ID(k).Type, 'malignancy', ~ID(k).IsBenign);
        end
      
    end
    
    if (options.showImages)       
        if (~plusRGBAnalysis) ; lines = [measured; estimated]; else;  lines = [measured; estimated; estRGB]; end   
        options.saveOptions.plotName = generateName('plot', options, ID(k));   
        plots('estimationComparison', 1, lines', 'Reflectance Estimation Comparison', 'Wavelength', wavelength, 'Method', ... 
            strcat(options.smoothingMatrixMethod, ' +  ', options.noiseType), ...
            'SaveOptions', options.saveOptions, 'LineNames', lineNames);
        
        pause(0.1)
    end
    
end
warning(w);

%% Export results
errors = GetErrorInfoStruct(rmse, nmse, pixelValueSelectionMethods, smoothingMatrixMethods, noiseTypes);

if (isNoComparison)
    % show NRMSE error per set
    
    nrmse = showNRMSE(rmse, Spectra, ID);
    plots('reconstructionPerformanceBars', 1, [], 'Reconstruction NRMSE', 'LineNames', {'benign', 'malignant'}, 'Performance', nrmse, 'SaveOptions', options.saveOptions)
    writeOutput(options, errFixingInfo, estimatedSpectra, hasRGBAnalysis && ~plusRGBAnalysis);
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
            if strcmp(options.dataset, 'saitama_v7_min_square_e')
                nms = {'spatial olympus'}; %{'sameForChannel 0.0001'};
            elseif strcmp(options.dataset, 'saitama_v7_min_region_e')
                 nms = {'spatial olympus'};
            elseif strcmp(options.dataset, 'saitama_v6_min_square_e')
                 nms = {'spatial olympus'};
            elseif strcmp(options.dataset, 'saitama_v6_min_region_e')
                 nms = {'spatial olympus'};
            else
                nms = {'spatial olympus'};
            end
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
            fprintf('NRMSE for malignancy %d  and state %s = %.4f\n',malignancy, state{1}, nrmseArray(i));               
        end
        nrmse = GetNRMSE(measured, rmse, strcmp({ID.Type}, state));
        fprintf('NRMSE for %s = %.4f\n', state{1}, nrmse);
    end
   nrmseArray = reshape(nrmseArray,[2, 3]); 
   
    for malignancy = 0:1
        nrmse = GetNRMSE(measured, rmse, [ID.IsBenign] ~= malignancy);
        fprintf('NRMSE overall for malignancy %d = %.4f\n',malignancy, nrmse);
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
    if exist(options.outName, 'file')
        if (hasRGBAnalysis)
            save(strrep(options.outName, '_rgb', ''), 'EstimatedRGBSpectra','-append');
        else
            save(strrep(options.outName, '_rgb', ''), 'EstimatedSpectra', '-append');
        end
    else
        if (hasRGBAnalysis)
            save(strrep(options.outName, '_rgb', ''), 'EstimatedRGBSpectra');
        else
            save(strrep(options.outName, '_rgb', ''), 'EstimatedSpectra');
        end
    end
end

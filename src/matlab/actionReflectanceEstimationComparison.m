%Start of reflectance estimation comparison
smmsFull = { 'Cor_All', 'Cor_Fixing', 'Cor_Sample', 'adaptive'}; %, 'Cor_SampleMalignancy', 'Cor_SampleMalignancyFixing', 'Cor_MalignancyFixing', 'Cor_SampleMalignancyFixing', 'markovian', 'Cor_Macbeth',};
pvsmsFull = {'green', 'rms', 'adjusted', 'extended', 'rgb'};
nmsFull = {'sameForChannel 0.0001', 'sameForChannel 0.001', 'sameForChannel 0.01',  'sameForChannel 0.1', 'sameForChannel 0.5',...
    'diffForChannel 0.0016, 0.0016, 0.0008, 0.0005, 0.0008, 0.0148, 0.0015', 'diffForChannel 0.16, 0.16, 0.08, 0.05, 0.08, 1.48, 0.15',...
    'SNR 25', 'SNR 20', 'SNR 15', 'SNR 5',  ...
    'spatial 0.002 0.04', 'spatial 0.005 0.5', 'spatial 0.002 0.01', 'fromOlympus'}; %, 'fromOlympus', 'spatial', 'spatial 0.015 0.015' 

%% Comparison settings
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
        nms = {'spatial olympus'}; %{'sameForChannel 0.0001'};
        plotType = '';
        
elseif contains(lower(options.action), 'extra') 
        pvsms = {'extended'};
        smms = {'Cor_All', 'adaptive'};
        nms = {'sameForChannel 0.0001', 'spatial 0.0004 0.002'};
        plotType = '';

else 
        error('Not implemented yet.')
end
isSingleMethod = contains(lower(options.action), 'preset') || contains(lower(options.action), 'simple');
hasRGBAnalysis = strcmp(options.pixelValueSelectionMethod, 'rgb') || contains(lower(options.action), 'rgb');

[lineNames, methodsN, pixelValueSelectionMethods, smoothingMatrixMethods, noiseTypes, smmsN, pvsmsN, nmsN] = getLineNames(pvsms, smms, nms);

%% Comparison
w = warning('off', 'all');

tic;
rmse = zeros(methodsN, msiN);
nmse = zeros(methodsN, msiN);
EstimatedSpectra = zeros(size(Spectra));
EstimatedRGBSpectra = zeros(size(Spectra));
nmseRGB = zeros(methodsN, msiN);
rmseRGB = zeros(methodsN, msiN);

errFixingInfo = struct('rmse', [], 'fixing', []);
errFixingInfoRGB = struct('rmse', [], 'fixing', []);
for k = 1:msiN
    % Give a name
    [options.saveOptions.plotName, sampleName] = generateName(options, 'plot', ID(k));   
    % Retrieve MSI data
    msi = MSIs{k};
    mask = Masks{k};
    measured = Spectra(k,:);
    if (hasRGB); rgb = readRGBImage(WhiteIs{k}); end
        
    % Estimation
    estimated = zeros(methodsN, length(wavelength));
    
    for j = 1:methodsN
        [options, optionsRgb]= updateOptions(options, j, pvsms, smms, nms);
        
        % In the case of 3channel RGB image
        if (hasRGB)
            [estRgb, rmseRGB(j, k), nmseRGB(j, k), ~] = reflectanceEstimation(rgb, mask, measured, ID(k), optionsRgb);
            if (isSingleMethod)
                EstimatedRGBSpectra(k,:) = estRgb;
                errFixingInfoRGB(k) = struct('rmse', rmseRGB(j, k), 'fixing', ID(k).Type);
            end
        end
        
        [estimated(j,:), rmse(j, k), nmse(j, k), ~] = reflectanceEstimation(msi, mask, measured, ID(k), options);                 
        if (isSingleMethod)
            EstimatedSpectra(k,:) = estimated(j,:);   
            errFixingInfo(k) = struct('rmse', rmse(j, k), 'fixing',  ID(k).Type);
        end   
    end
    
    if (options.showImages)       
        if (~hasRGB) ; lines = [measured; estimated]; else;  lines = [measured; estimated; estRgb]; end;        
        plots('estimationComparison', 1, lines', sampleName, 'Wavelength', wavelength, 'Method', ... 
            strcat(options.smoothingMatrixMethod, ' +  ', options.noiseType), ...
            'SaveOptions', options.saveOptions, 'LineNames', lineNames);
        
        pause(0.1)
    end
    
end
warning(w);

%% Export results
errors = GetErrorInfoStruct(rmse, nmse, pixelValueSelectionMethods, smoothingMatrixMethods, noiseTypes)
if (hasRGB);  errorsRGB = GetErrorInfoStruct(rmseRGB, nmseRGB, pixelValueSelectionMethods, smoothingMatrixMethods, noiseTypes); end

fprintf('Minimum rmse = %.5f\n', min(mean(rmse, 2)));

if (isSingleMethod)
    writetable(struct2table(errFixingInfo), strcat('../../logs/reflectance_estimation_comparison_', options.dataset,'_log.xlsx'), 'Sheet', 'MSI');
    writetable(struct2table(errFixingInfoRGB), strcat('../../logs/reflectance_estimation_comparison_', options.dataset,'_log.xlsx'), 'Sheet', 'RGB');
    if exist(options.outName, 'file')
        save(strrep(options.outName, '_rgb', ''), 'EstimatedSpectra', 'EstimatedRGBSpectra', 'rmse', '-append');
    else
        save(strrep(options.outName, '_rgb', ''), 'EstimatedSpectra', 'EstimatedRGBSpectra', 'rmse');
    end
    
else
    writetable(struct2table(errors), strcat('../../logs/reflectance_estimation_comparison_', options.dataset,'_log.xlsx'));

    options.saveOptions.plotName = generateName(options, ['ComparisonRmse', plotType]);
    plots('methodErrors', 3, [], ['Rmse', plotType], 'Errors', errors, 'SaveOptions', options.saveOptions)

    options.saveOptions.plotName = generateName(options, ['ComparisonNmse', plotType]);
    plots('methodErrors', 4, [], ['Nmse', plotType], 'Errors', errors, 'SaveOptions', options.saveOptions)
    
    save(generateName(options, ['ComparisonErrors', plotType, '.mat']), 'errors');
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

function [updatedOptions] = updateOptions(options, j, pvsms, smms, nms)
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

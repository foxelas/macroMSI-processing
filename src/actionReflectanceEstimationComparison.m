%Start of reflectance estimation comparison
smmsFull = { 'Cor_All', ... 'Cor_Malignancy', 'Cor_Fixing', 'Cor_Sample', ...
     'adaptive'}; %, 'Cor_SampleMalignancy', 'Cor_SampleMalignancyFixing', 'Cor_MalignancyFixing', 'Cor_SampleMalignancyFixing', 'markovian', 'Cor_Macbeth',};
pvsmsFull = {'green', 'rms', 'adjusted', 'extended'};
nmsFull = {'sameForChannel 0.001', 'diffForChannel [0.0031, 0.0033, 0.0030, 0.0031, 0.0032, 0.0029, 0.0024]', ...
    'givenSNR 17.5dB', 'spatial 0.02 0.04'}; %, 'fromOlympus', 'spatial', 'spatial 0.015 0.015' 


%% Comparison settings
if contains(lower(options.action), 'matrixsystem')
        smms = smmsFull;
        pvsms = pvsmsFull;
        plotType = 'MatrixSystemAvg';        
        nms = {'givenSNR'}; 
        
elseif contains(lower(options.action), 'matrixnoise')
        smms = smmsFull;
        pvsms = {'extended'};
        plotType = 'MatrixNoiseAvg';        
        nms = nmsFull; 
        
elseif contains(lower(options.action), 'smoothing matrix') || contains(lower(options.action), 'matrix')
        smms = smmsFull;
        pvsms = {'extended'};
        plotType = 'MatrixAvgMinMax';
        nms = {'givenSNR'};
        
elseif contains(lower(options.action), 'color system') || contains(lower(options.action), 'system')   
        smms = {'Cor_Malignancy'};
        pvsms = pvsmsFull;
        plotType = 'SystemAvgMinMax';
        nms = {'givenSNR'};

elseif contains(lower(options.action), 'noise model') || contains(lower(options.action), 'noise')   
        pvsms = {'extended'};
        smms = {'Cor_SampleMalignancy'};
        plotType = 'NoiseAvgMinMax'; %Or 'NoiseAvgMinMaxWithNone'
        nms = nmsFull;

elseif contains(lower(options.action), 'simple') 
        pvsms = {'extended'};
        smms =  {options.smoothingMatrixMethod};
        nms = {options.noiseType};
        estimatedSpectrumStruct = struct('Name', {}, 'Index', [], 'Spectrum', []);
        
elseif contains(lower(options.action), 'preset')
        pvsms = {'extended'};
        smms =  {'Cor_SampleMalignancy'}; %{'Cor_Malignancy'}; 
        nms = {'givenSNR'};
        estimatedSpectrumStruct = struct('Name', {}, 'Index', [], 'Spectrum', []);    
     
else 
        error('Not implemented yet.')
end

smmsN = numel(smms);
pvsmsN = numel(pvsms);
nmsN = length(nms);
methodsN = smmsN * pvsmsN * nmsN;

optionsSelection = repmat(options, methodsN, 1);
for l = 1:pvsmsN
    for n = 1:smmsN
        for m = 1:nmsN
            optionsSelection((m - 1) * pvsmsN * smmsN + (n - 1) * pvsmsN + l).pixelValueSelectionMethod = pvsms{l};
            optionsSelection((m - 1) * pvsmsN * smmsN + (n - 1) * pvsmsN + l).smoothingMatrixMethod = smms{n};
            optionsSelection((m - 1) * pvsmsN * smmsN + (n - 1) * pvsmsN + l).noiseType = nms{m};
        end
    end
end
rmse = zeros(methodsN, msiN);
nmse = zeros(methodsN, msiN);

%% Comparison
% range = 146:187; %; %1:msiN; [66, 165, 175, 236, 297]
for k = 1:msiN
    % Give a name
    [options.saveOptions.plotName, sampleName] = generateName(options, 'plot', ID(k));
    
    % Retrieve MSI data
    g = MSIStruct(k);
    
    % Retrieve measured spectrum
    measured = interp1(380:780, measuredSpectrumStruct(k).Spectrum, wavelength, 'nearest');
    
    %% Estimation
    est = zeros(length(wavelength), methodsN);
    lineNames = {'Center \lambda', 'Measured'};
    for l = 1:pvsmsN
        for n = 1:smmsN
            for m = 1:nmsN
            
                if strcmp(optionsSelection(l).pixelValueSelectionMethod, 'rgb')
                    tempRGB = whiteStruct(k).MSI;
                    g.MSI = reshape(tempRGB, [3, size(tempRGB, 1), size(tempRGB, 2)]);
                end

                j = (m - 1) * pvsmsN * smmsN + (n - 1) * pvsmsN + l ;

                % optionsSelection(j).coeff = squeeze(coeff(k, min(find(strcmp(pvsms , optionsSelection(j).pixelValueSelectionMethod)),3), :))';
                [est(:, j), rmse(j, k), nmse(j, k)] = reflectanceEstimation(g, measured, ID(k), optionsSelection(j));
                if (pvsmsN < 2)
                    lineNames{j+2} = strcat('Est-', smms{n});
                elseif (smmsN < 2)
                    lineNames{j+2} = strcat('Est-', pvsms{l});
                else 
                    lineNames{j+2} = strcat('Est-', pvsms{l}, '|', smms{n});
                end
                if (nmsN > 1)
                    lineNames{j+2} = strcat('Est-', nms{m});
                end

                if contains(lower(options.action), 'preset') || contains(lower(options.action), 'simple')
                    estimatedSpectrumStruct(k) = struct('Name', sampleName, 'Index', MSIStruct(k).Index, 'Spectrum', est(:, j));
                end
                
            end
        end
    end
    
    if (options.showImages)
        plots('estimationComparison', 2, [measured, est], sampleName, 'wavelength', wavelength, 'method', options.pixelValueSelectionMethod, ...
            'saveOptions', options.saveOptions, 'lineNames', lineNames);
        pause(0.1)
    end
    
end

%% Export results
rmse = rmse(:, rmse(1, :) ~= 0);
nmse = nmse(:, nmse(1, :) ~= 0);
errorData = [mean(rmse, 2), max(rmse, [], 2), min(rmse, [], 2), mean(nmse, 2), max(nmse, [], 2), min(nmse, [], 2)];
errors = struct('avgrmse', mean(rmse, 2), 'minrmse', min(rmse, [], 2), 'maxrmse', max(rmse, [], 2), 'stdrmse', std(rmse, [], 2), 'avgnmse', mean(nmse, 2), 'minnmse', min(nmse, [], 2), 'maxnmse', max(nmse, [], 2), 'stdnmse', std(nmse, [], 2), 'options', optionsSelection);

minError = min(mean(rmse, 2));
fprintf('Minimum rmse = %.5f\n', minError);

if contains(lower(options.action), 'preset') || contains(lower(options.action), 'simple')
    out.EstimatedSpectrumStruct = estimatedSpectrumStruct;
else
    options.saveOptions.plotName = generateName(options, ['ComparisonRmse', plotType]);
    plots('methodErrors', 3, [], ['Rmse', plotType], 'errors', errors, 'saveOptions', options.saveOptions)

    options.saveOptions.plotName = generateName(options, ['ComparisonNmse', plotType]);
    plots('methodErrors', 4, [], ['Nmse', plotType], 'errors', errors, 'saveOptions', options.saveOptions)

    errorfilename = generateName(options, ['ComparisonErrors', plotType, '.mat']);
    save(errorfilename, 'errors');
end

%End of reflectance estimation comparison

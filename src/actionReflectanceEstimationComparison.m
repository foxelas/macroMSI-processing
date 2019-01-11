%Start of reflectance estimation comparison
smmsFull = {'markovian', 'Cor_Macbeth', 'Cor_All', 'Cor_Malignancy',...
    'Cor_Fixing', 'Cor_MalignancyFixing', 'Cor_Sample'}; %, 'Cor_SampleMalignancyFixing'};
pvsmsFull = {'green', 'rms', 'adjusted', 'extended'};
nmsFull = {'none', 'white gaussian 10^{-5}', 'independent 10^{-3}', 'independent 10^{-5}', 'givenSNR 15dB', 'givenSNR 25dB', 'givenSNR 30dB', 'givenSNR 40dB', 'fromOlympus', 'spatial'};

%% Comparison settings
if contains(action, 'matrixsystem')
        smms = smmsFull;
        pvsms = pvsmsFull;
        plotType = 'MatrixSystemAvg';        
        nms = {'givenSNR'}; 
        
elseif contains(action, 'matrixnoise')
        smms = smmsFull;
        pvsms = {'extended'};
        plotType = 'MatrixNoiseAvg';        
        nms = nmsFull; 
        
elseif contains(action, 'smoothing matrix') || contains(action, 'matrix')
        smms = smmsFull;
        pvsms = {'extended'};
        plotType = 'MatrixAvgMinMax';
        nms = {'givenSNR'};
        
elseif contains(action, 'color system') || contains(action, 'system')   
        smms = {'Cor_Malignancy'};
        pvsms = pvsmsFull;
        plotType = 'SystemAvgMinMax';
        nms = {'givenSNR'};

elseif contains(action, 'noise model') || contains(action, 'noise')   
        pvsms = {'extended'};
        smms = {'Cor_Malignancy'};
        plotType = 'NoiseAvgMinMax'; %Or 'NoiseAvgMinMaxWithNone'
        nms = nmsFull;

elseif contains(action, 'simple') 
        pvsms = {'extended'};
        smms = {'Cor_Malignancy'};
        nms = {options.noiseType}; %{'givenSNR'};
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
load(fullfile(options.systemdir, 'patches.mat'));
%% Comparison
% range = 146:187; %; %1:msiN; [66, 165, 175, 236, 297]
for k = 1:5 %1:msiN
    % Give a name
    [options.saveOptions.plotName, sampleName] = generateName(options, 'plot+save', data(ID(k).Representative), ID(k));
    
    % Retrieve MSI data
%     g = MSIStruct(k).MSI;
    g = patches(k);

    
    % Retrieve measured spectrum
    measured = interp1(380:780, measuredSpectrumStruct(k).Spectrum, wavelength, 'nearest');
    
    %% Estimation
    est = zeros(length(wavelength), methodsN);
    lineNames = {'Center \lambda', 'Measured'};
    for l = 1:pvsmsN
        for n = 1:smmsN
            for m = 1:nmsN
            
                if strcmp(optionsSelection(l).pixelValueSelectionMethod, 'rgb')
                    g = whiteStruct(k).MSI;
                    g = reshape(g, [3, size(g, 1), size(g, 2)]);
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

                if contains(action, 'simple') 
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
% fprintf('Minimum rmse = %.5f\n', minError);

if contains(action, 'simple') 
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

tic;

%Start of reflectance estimation comparison

[optionsSelection, lineNames, plotType] = getComparisons(options);
methodsN = length(optionsSelection);
estimatedSpectrumStruct = struct('Name', {}, 'Index', [], 'Spectrum', []);
rmse = zeros(methodsN, msiN);
nmse = zeros(methodsN, msiN);
newCoordinates = zeros(msiN, 2);

%% Comparison
measured = cellfun(@(x) interp1(380:780, x, wavelength, 'nearest'), {measuredSpectrumStruct.Spectrum}, 'un', 0);
for k = 1:msiN
    % Give a name
    [options.saveOptions.plotName, sampleName] = generateName(options, 'plot', ID(k));
    
    % Retrieve MSI data
    g = MSIStruct(k);
    
    if isfield(options,'single')
        [r, c] = find(g.MaskI);
        x = ID(k).Originx - min(c) + 1;
        y = ID(k).Originy - min(r) + 1;
        g.MSI = g.MSI(:,y, x,:);
        g.Mask = 1;
    end
        
    %% Estimation
    est = zeros(length(wavelength), methodsN);
    for j = 1:methodsN
            
        if strcmp(optionsSelection(j).pixelValueSelectionMethod, 'rgb')
            tempRGB = whiteStruct(k).MSI;
            g.MSI = reshape(tempRGB, [3, size(tempRGB, 1), size(tempRGB, 2)]);
        end

        [est(:, j), rmse(j, k), nmse(j, k), minIdx] = reflectanceEstimation(g, measured{k}, ID(k), optionsSelection(j));

        if contains(lower(options.action), 'preset') || contains(lower(options.action), 'simple')
            
            estimatedSpectrumStruct(k) = struct('Name', sampleName, 'Index', MSIStruct(k).Index, 'Spectrum', est(:, j));
            
            [r, c] = find(g.MaskI);
            newCoordinates(k, 1:2) = [min(c) + minIdx(2) - 1, min(r)+ minIdx(1) - 1];
            if (options.showImages)
                seg = readMSI({data(ID(k).Data).File});
                coordinates = [newCoordinates(k, 1:2); ID(k).Originx, ID(k).Originy];
                plots('cropped', 2, 'Image', seg.whiteReference + cat(3, 0.1*g.MaskI, 0.05*g.MaskI, 0.1*g.MaskI), 'Coordinates', coordinates );
            end            
        end    
        
    end
    
    if (options.showImages)
        plots('estimationComparison', 1, [measured{k}, est], sampleName, 'Wavelength', wavelength, 'Method', ...
            strcat(optionsSelection(j).smoothingMatrixMethod, ' +  ', optionsSelection(j).noiseType), ...
            'SaveOptions', options.saveOptions, 'LineNames', lineNames);
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
    out.newCoordinates = newCoordinates;
    
else
    options.saveOptions.plotName = generateName(options, ['ComparisonRmse', plotType]);
    plots('methodErrors', 3, [], ['Rmse', plotType], 'Errors', errors, 'SaveOptions', options.saveOptions)

    options.saveOptions.plotName = generateName(options, ['ComparisonNmse', plotType]);
    plots('methodErrors', 4, [], ['Nmse', plotType], 'Errors', errors, 'SaveOptions', options.saveOptions)

    errorfilename = generateName(options, ['ComparisonErrors', plotType, '.mat']);
    save(errorfilename, 'errors');
end
t = toc
%End of reflectance estimation comparison

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [optionsSelection, lineNames, plotType] = getComparisons(options)

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
            plotType = '';

    elseif contains(lower(options.action), 'preset')
            pvsms = {'extended'};
            smms =  {'adaptive'}; %{'Cor_Malignancy'}; 
            nms = {'givenSNR'};
            plotType = '';
            
    else 
            error('Not implemented yet.')
    end

    smmsN = numel(smms);
    pvsmsN = numel(pvsms);
    nmsN = length(nms);
    methodsN = smmsN * pvsmsN * nmsN;

    lineNames = {'Center \lambda', 'Measured'};
    optionsSelection = repmat(options, methodsN, 1);
    for l = 1:pvsmsN
        for n = 1:smmsN
            for m = 1:nmsN
                optionsSelection((m - 1) * pvsmsN * smmsN + (n - 1) * pvsmsN + l).pixelValueSelectionMethod = pvsms{l};
                optionsSelection((m - 1) * pvsmsN * smmsN + (n - 1) * pvsmsN + l).smoothingMatrixMethod = smms{n};
                optionsSelection((m - 1) * pvsmsN * smmsN + (n - 1) * pvsmsN + l).noiseType = nms{m};
                
                j = (m - 1) * pvsmsN * smmsN + (n - 1) * pvsmsN + l ;    
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
            end
        end
    end

end
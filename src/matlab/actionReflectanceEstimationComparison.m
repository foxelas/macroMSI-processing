%Start of reflectance estimation comparison

smmsFull = { 'Cor_All', 'Cor_Fixing', 'Cor_Sample', 'adaptive'}; %, 'Cor_SampleMalignancy', 'Cor_SampleMalignancyFixing', 'Cor_MalignancyFixing', 'Cor_SampleMalignancyFixing', 'markovian', 'Cor_Macbeth',};
pvsmsFull = {'green', 'rms', 'adjusted', 'extended', 'rgb'};
nmsFull = {'sameForChannel 0.0001', 'diffForChannel 0.0016, 0.0016, 0.0008, 0.0005, 0.0008, 0.0148, 0.0015', 'SNR 25', 'spatial 0.002 0.04'}; %, 'fromOlympus', 'spatial', 'spatial 0.015 0.015' 

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
        nms = {'spatial 0.4 0.2'}; %{'sameForChannel 0.0001'};
        plotType = '';
        
elseif contains(lower(options.action), 'extra') 
        pvsms = {'extended'};
        smms = {'Cor_All', 'adaptive'};
        nms = {'sameForChannel 0.0001', 'spatial 0.0004 0.002'};
        plotType = '';

else 
        error('Not implemented yet.')
end

smmsN = numel(smms);
pvsmsN = numel(pvsms);
nmsN = length(nms);
methodsN = smmsN * pvsmsN * nmsN;

lineNames = {'Channel \lambda', 'Measured'};
for l = 1:pvsmsN
    for n = 1:smmsN
        for m = 1:nmsN
            j = sub2ind([ pvsmsN , smmsN, nmsN],  l, n, m);     
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

if contains(lower(options.action), 'rgb')
    lineNames{end + 1} = 'RGB Estimate';
end

rmse = zeros(methodsN, msiN);
nmse = zeros(methodsN, msiN);
newCoordinates = zeros(msiN, 2);
isSingleMethod = contains(lower(options.action), 'preset') || contains(lower(options.action), 'simple');

%% Comparison
w = warning('off', 'all');

tic;
EstimatedSpectra = zeros(size(Spectra));
for k = 1:msiN
    % Give a name
    [options.saveOptions.plotName, sampleName] = generateName(options, 'plot', ID(k));   
    % Retrieve MSI data
    msi = MSIs{k};
    mask = Masks{k};
    measured = Spectra(k,:);
        
    % Estimation
    estimated = zeros(methodsN, length(wavelength));
    
    for j = 1:methodsN
        [l,n,m] = ind2sub([ pvsmsN , smmsN, nmsN], j);
        options.pixelValueSelectionMethod = pvsms{l};
        options.smoothingMatrixMethod = smms{n};
        options.noiseType = nms{m};   
        
        % In the case of 3channel RGB image
        if strcmp(options.pixelValueSelectionMethod, 'rgb')
            tempRGB = WhiteIs{k};
            msi = reshape(tempRGB, [3, size(tempRGB, 1), size(tempRGB, 2)]);
        end

%         if isfield(options,'SinglePixel')
%             maskI = MaskIs{k};
%             [r, c] = find(maskI);
%             x = ID(k).Originx - min(c) + 1;
%             y = ID(k).Originy - min(r) + 1;
%             msi = msi(:,y, x,:);
%             mask = 1;
%         end
        
        [estimated(j,:), rmse(j, k), nmse(j, k), minIdx] = reflectanceEstimation(msi, mask, measured, ID(k), options);
        
        if (isSingleMethod)
            EstimatedSpectra(k,:) = estimated(j,:); 
%             maskI = MaskIs{k};
%             [r, c] = find(maskI);
%             newCoordinates(k, 1:2) = [min(c) + minIdx(2) - 1, min(r)+ minIdx(1) - 1];          
        end    
        
    end
    
    if (options.showImages)
        
        lines = [measured; estimated];
        
        % For comparison with RGB estimation
        if contains(lower(options.action), 'rgb')
            optionsRgb = options;
            optionsRgb.pixelValueSelectionMethod = 'rgb';
            tempRGB = WhiteIs{k};
            rgb = reshape(tempRGB, [3, size(tempRGB, 1), size(tempRGB, 2)]);
            estRgb = reflectanceEstimation(rgb, mask, measured, ID(k), optionsRgb)'; 
            lines = [measured; estimated; estRgb];
        end
            
        plots('estimationComparison', 1, lines', sampleName, 'Wavelength', wavelength, 'Method', ... 
            strcat(options.smoothingMatrixMethod, ' +  ', options.noiseType), ...
            'SaveOptions', options.saveOptions, 'LineNames', lineNames);
        
%         [r, c] = find(maskI);
%         hold on 
%         scatter( [450, 465, 505, 525, 575, 605, 630], squeeze( raw2msi(msi(:, minIdx(1), minIdx(2), :), 'adjusted')) * 25, 'mo');
%         hold off
        if false %(isSingleMethod)
            [~, whiteReference] = readMSI({data([data.MsiID] == ID(k).MsiID).File});
            coordinates = [newCoordinates(k, 1:2); ID(k).Originx, ID(k).Originy];
            plots('cropped', 2, 'Image', whiteReference + cat(3, 0.1*maskI, 0.05*maskI, 0.1*maskI), 'Coordinates', coordinates );
        end
        
        pause(0.1)
    end
    
end
warning(w);

%% Export results
% rmse = rmse(:, rmse(1, :) ~= 0);
% nmse = nmse(:, nmse(1, :) ~= 0);
errorData = [mean(rmse, 2), max(rmse, [], 2), min(rmse, [], 2), mean(nmse, 2), max(nmse, [], 2), min(nmse, [], 2)];

pixelValueSelectionMethods = cell(methodsN,1);
smoothingMatrixMethods = cell(methodsN,1);
noiseTypes = cell(methodsN,1);

for j = 1:methodsN
    [l,n,m] = ind2sub([ pvsmsN , smmsN, nmsN], j);
    pixelValueSelectionMethods{j} = pvsms{l};
    smoothingMatrixMethods{j} = smms{n};
    noiseTypes{j} = nms{m};   
end
errors = struct('avgrmse', num2cell(mean(rmse, 2)), 'minrmse', num2cell(min(rmse, [], 2)), 'maxrmse', num2cell(max(rmse, [], 2)), ...
                'stdrmse', num2cell(std(rmse, [], 2)), 'avgnmse', num2cell(mean(nmse, 2)), 'minnmse', num2cell(min(nmse, [], 2)), ...
                'maxnmse', num2cell(max(nmse, [], 2)), 'stdnmse', num2cell(std(nmse, [], 2)), 'pixelValueSelectionMethod', pixelValueSelectionMethods, ...
                'smoothingMatrixMethod', smoothingMatrixMethods, 'noiseType', noiseTypes);

minError = min(mean(rmse, 2));
fprintf('Minimum rmse = %.5f\n', minError);

if contains(lower(options.action), 'preset') || contains(lower(options.action), 'simple')
    fprintf('RMSE: avg=%.5f, std=%.5f, min=%.5f, max=%.5f\n',  mean(rmse), std(rmse), min(rmse), max(rmse));
    fileID = fopen('..\..\logs\preset_log.txt','a+');
    fprintf(fileID,'%s \n', options.dataset);
    fprintf(fileID,'RMSE: avg=%.5f, std=%.5f, min=%.5f, max=%.5f\n',  mean(rmse), std(rmse), min(rmse), max(rmse));
    fclose(fileID);

    goodSampleIndexes_a = arrayfun(@(x) ~any(x == [233, 274, 310, 349, 378, 512]), (1:length(rmse))');
    goodSampleIndexes_b = (rmse < mean(rmse))';
    goodSampleIndexes = goodSampleIndexes_a & goodSampleIndexes_b;
    if  ~isfile(options.outName)
        dateCreated = datetime();
        save(options.outName, 'dateCreated');
    end
    save(options.outName, 'EstimatedSpectra', 'rmse', 'goodSampleIndexes', '-append');
%     out.NewCoordinates = newCoordinates;
    
else
    options.saveOptions.plotName = generateName(options, ['ComparisonRmse', plotType]);
    plots('methodErrors', 3, [], ['Rmse', plotType], 'Errors', errors, 'SaveOptions', options.saveOptions)

    options.saveOptions.plotName = generateName(options, ['ComparisonNmse', plotType]);
    plots('methodErrors', 4, [], ['Nmse', plotType], 'Errors', errors, 'SaveOptions', options.saveOptions)

    errorfilename = generateName(options, ['ComparisonErrors', plotType, '.mat']);
    save(errorfilename, 'errors');
end
t2 = toc; 
fprintf('Run time time %.5f\n', t2);

%End of reflectance estimation comparison


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


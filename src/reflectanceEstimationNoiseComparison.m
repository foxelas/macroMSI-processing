%Start of reflectance estimation noise model comparison

nms = {'none', 'white gaussian 10^{-5}', 'independent 10^{-3}', 'independent 10^{-5}', 'givenSNR 15dB', 'givenSNR 25dB', 'givenSNR 30dB', 'givenSNR 40dB','fromOlympus'}; 
nmsN = length(nms);
options.pixelValueSelectionMethod = 'extended';
options.smoothingMatrixMethod = 'KCor all specimen';
optionsSelection = repmat(options, nmsN, 1);
for n = 1:nmsN
    optionsSelection(n).noiseType = nms{n};
end
optionsTmp = optionsSelection;
for n = 1:nmsN
    if contains(optionsSelection(n).noiseType, 'independent')
        attr = strsplit(optionsSelection(n).noiseType, {' ', '^{', '}'});
        optionsTmp(n).noiseOrder = str2double(attr{2})^str2double(attr{3});
        optionsTmp(n).noiseType = 'independent';
    end
    if contains(optionsSelection(n).noiseType, 'givenSNR')
        attr = strsplit(optionsSelection(n).noiseType, {' ', 'dB'});
        optionsTmp(n).snr = str2double(attr{2});
        optionsTmp(n).noiseType = 'givenSNR';
    end
    if contains(optionsSelection(n).noiseType, 'white gaussian')
        attr = strsplit(optionsSelection(n).noiseType, {' ', '^{', '}'});
        optionsTmp(n).noiseOrder = str2double(attr{3})^str2double(attr{4});
        optionsTmp(n).noiseType = 'white gaussian';
    end
end

rmse = zeros(nmsN, msiN); 
nmse = zeros(nmsN, msiN); 

% range = 146:187; %; %1:msiN; [66, 165, 175, 236, 297]
for k = 1:msiN
    % Give a name
    [options.saveOptions.plotName, sampleName] = generateName(options, 'plot+save',  data(ID(k).Representative), ID(k));    

    % Retrieve MSI data 
     g = MSIStruct(k).MSI; 

    % Retrieve measured spectrum
    measured = interp1(380:780, measuredSpectrumStruct(k).Spectrum, wavelength, 'nearest');   
    
   %% Estimation
     est = zeros(length(wavelength), nmsN);
     lineNames = {'MS center \lambda', 'Measured'};     
%      coeffs = squeeze(coeff(k, 3, :))'; % because of 'extended'
    
    for n = 1:nmsN
        % Retrieve correction coefficients
        optionsTmp(n).coeff = coeffs;     
        
        [est(:, n), rmse(n,k), nmse(n,k)] = reflectanceEstimation( g, measured, ID(k), optionsTmp(n)); 
        lineNames{n+2} = strcat('Est-', nms{n});
    end
   
     if (options.showImages)
        plots( 'estimationComparison', 2, [measured, est], sampleName, 'wavelength', wavelength, 'method', options.pixelValueSelectionMethod, ...
            'saveOptions', options.saveOptions, 'lineNames', lineNames);
     end
end

%% Export errors 
rmse = rmse(:,rmse(1,:) ~= 0);
nmse = nmse(:,nmse(1,:) ~= 0);
errorData = [mean(rmse, 2), max(rmse,[], 2), min(rmse, [], 2), mean(nmse, 2), max(nmse,[], 2), min(nmse, [], 2)];
errors = struct('avgrmse', mean(rmse, 2), 'minrmse', min(rmse,[], 2), 'maxrmse',max(rmse, [], 2), 'stdrmse', std(rmse,[], 2), 'avgnmse', mean(nmse, 2), 'minnmse', min(nmse, [], 2), 'maxnmse', max(nmse,[], 2), 'stdnmse', std(nmse,[], 2), 'options', optionsSelection );

plots('methodErrors', 1, [], 'RmseNoiseAvgMinMaxWithNone', 'errors', errors);
plots('methodErrors', 2, [], 'NmseNoiseAvgMinMaxWithNone', 'errors', errors);
save(strcat(options.saveOptions.savedir, 'errorsNoiseAvgMinMax.mat') , 'errors');

%without 'none' option
errors = structfun(@(x) x(3:end),errors, 'UniformOutput', false);

options.saveOptions.plotName = generateName(options, 'override', [], [], 'ComparisonRmseNoiseAvgMinMax');
plots('methodErrors', 3, [], 'RmseNoiseAvgMinMax', 'errors', errors, 'saveOptions', options.saveOptions)

options.saveOptions.plotName = generateName(options, 'override', [], [], 'ComparisonNmseNoiseAvgMinMax');
plots('methodErrors', 4, [], 'NmseNoiseAvgMinMax', 'errors', errors, 'saveOptions', options.saveOptions)

%End of reflectance estimation noise model comparison
        
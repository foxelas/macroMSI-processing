%Start of reflectance estimation system comparison

pvsms = pixelValueSelectionMethods; 
methodsN = length(pvsms);
options.smoothingMatrixMethod = 'KCor all same malignancy';
optionsSelection = repmat(options, methodsN, 1);
for l = 1:methodsN
        optionsSelection(l).pixelValueSelectionMethod = pvsms{l};
end
rmse = zeros(methodsN, msiN); 
nmse = zeros(methodsN, msiN); 

% range = 146:187; %1:msiN; % [57, 125, 215, 278] % [66, 165, 175, 236, 297]
for k = 1:msiN
    % Give a name
    [options.saveOptions.plotName, sampleName] = generateName(options, 'plot+save',  data(ID(k).Representative), ID(k));    

    % Retrieve MSI data 
     g = MSIStruct(k).MSI; 

    % Retrieve measured spectrum
    measured = interp1(380:780, measuredSpectrumStruct(k).Spectrum, wavelength, 'nearest');   
    
     %% Estimation
     est = zeros(length(wavelength), methodsN);
     for l = 1:methodsN
         
        if strcmp(optionsSelection(l).pixelValueSelectionMethod, 'rgb')
            g = whiteStruct(k).MSI; 
            g = reshape(g, [3, size(g,1), size(g,2)]);
        end
        
%         optionsSelection(l).coeff = squeeze(coeff(k, min(find(strcmp(pixelValueSelectionMethods , optionsSelection(l).pixelValueSelectionMethod)),3), :))';
        
        [est(:, l), rmse(l,k), nmse(l,k)] = reflectanceEstimation( g, measured, ID(k), optionsSelection(l)); 
     end
     if (options.showImages)
             plots( 'estimationComparison', 2, [measured, est], sampleName, 'wavelength', wavelength, 'method', options.smoothingMatrixMethod, ...
                 'saveOptions', options.saveOptions); %, 'MSIreflectances', msireflectances);  
     end
end

%% Export errors 
rmse = rmse(:,rmse(1,:) ~= 0);
nmse = nmse(:,nmse(1,:) ~= 0);
errors = struct('avgrmse', mean(rmse, 2), 'minrmse', min(rmse,[], 2), 'maxrmse',max(rmse, [], 2), 'stdrmse', std(rmse,[], 2), 'avgnmse', mean(nmse, 2), 'minnmse', min(nmse, [], 2), 'maxnmse', max(nmse,[], 2), 'stdnmse', std(nmse,[], 2), 'options', optionsSelection );

options.saveOptions.plotName = generateName(options, 'override', [], [], 'ComparisonRmseSystemAvgMinMax');
plots('methodErrors', 3, [], 'RmseSystemAvgMinMax', 'errors', errors, 'saveOptions', options.saveOptions)

options.saveOptions.plotName = generateName(options, 'override', [], [], 'ComparisonNmseSystemAvgMinMax');
plots('methodErrors', 4, [], 'NmseSystemAvgMinMax', 'errors', errors, 'saveOptions', options.saveOptions)

errorfilename = generateName(options, 'override', [], [], 'ComparisonErrorsSystemAvgMinMax.mat');
save(errorfilename, 'errors');

% end of reflectance estimation system comparison

        
%Start of reflectance estimation matrix system comparison

smms = {'Markovian 0.99', 'KCor all data', 'KCor macbeth', 'KCor all specimen', ...
          'KCor all same malignancy', 'KCor all same fixation' 'KCor same malignancy', 'KCor same malignancy, fixation'}; 
smmsN = length(smms);
pvsms = {'green',    'rms',    'adjusted',    'extended'}; 
pvsmsN = length(pvsms);
methodsN = smmsN * pvsmsN;

optionsSelection = repmat(options, methodsN, 1);
for l = 1:pvsmsN
    for n = 1:smmsN
        optionsSelection((l-1)*smmsN+n).pixelValueSelectionMethod = pvsms{l};
        optionsSelection((l-1)*smmsN+n).smoothingMatrixMethod = smms{n};
    end
end
rmse = zeros(methodsN, msiN); 
nmse = zeros(methodsN, msiN); 

% range = 146:187; %; %1:msiN; [66, 165, 175, 236, 297]
for k = 1:msiN
    % Give a name
    [options.saveOptions.plotName, sampleName] = generateName(options, 'plot+save',  data(ID(k).Representative), ID(k));    

    % Retrieve MSI data 
     g = MSIStruct(k).MSI;

    % Retrieve measured spectrum
    measured = interp1(380:780, measuredSpectrumStruct(k).Spectrum, wavelength, 'nearest'); 
    
   %% Estimation
     est = zeros(length(wavelength), methodsN);
     lineNames = {'MS center \lambda', 'Measured'};
     for l = 1:pvsmsN        
        for n = 1:smmsN   
            % Retrieve correction coefficients
            j = (l-1)*smmsN+n;
%             optionsSelection(j).coeff = squeeze(coeff(k, min(find(strcmp(pvsms , optionsSelection(j).pixelValueSelectionMethod)),3), :))';     
            [est(:, j), rmse(j,k), nmse(j,k)] = reflectanceEstimation( g, measured, ID(k), optionsSelection(j)); 
            lineNames{j+2} = strcat('Est-', pvsms{l}, '|', smms{n});
         end
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

options.saveOptions.plotName = generateName(options, 'override', [], [], 'ComparisonRmseMatrixSystemAvg');
plots('methodErrors', 3, [], 'RmseMatrixSystemAvg', 'errors', errors, 'saveOptions', options.saveOptions)

options.saveOptions.plotName = generateName(options, 'override', [], [], 'ComparisonNmseMatrixSystemAvg');
plots('methodErrors', 4, [], 'NmseMatrixSystemAvg', 'errors', errors, 'saveOptions', options.saveOptions)

errorfilename = generateName(options, 'override', [], [], 'ComparisonErrorsMatrixSystemAvg.mat');
save(errorfilename, 'errors');

%End of reflectance estimation system comparison

        
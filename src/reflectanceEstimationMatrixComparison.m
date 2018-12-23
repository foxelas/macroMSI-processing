%Start of reflectance estimation matrix comparison

smms = {'Markovian 0.99', 'KCor all data', 'KCor macbeth', 'KCor all specimen', ...
          'KCor all same malignancy', 'KCor all same fixation' 'KCor same malignancy', 'KCor same malignancy, fixation'}; 

methodsN = length(smms);

optionsSelection = repmat(options, methodsN, 1);
for l = 1:methodsN
    optionsSelection(l).smoothingMatrixMethod = smms{l};
end
rmse = zeros(methodsN, msiN); 
nmse = zeros(methodsN, msiN); 
% range = 146:187; %1:msiN; [66, 165, 175, 236, 297]
for k = 1:msiN
    % Give a name
    [options.saveOptions.plotName, sampleName] = generateName(options, 'plot+save',  data(ID(k).Representative), ID(k));    

    % Retrieve MSI data 
     g = MSIStruct(k).MSI; % eval(strcat('in.g_',  sampleName)); 

    % Retrieve measured spectrum
    measured = interp1(380:780, measuredSpectrumStruct(k).Spectrum, wavelength, 'nearest');   %interp1(380:780, eval(strcat('in.s_',  generateName([], 'csvfile', [], ID(k)) )), wavelength, 'nearest');   

    %% Estimation
      est = zeros(length(wavelength), methodsN);
     for l = 1:methodsN
         % Retrieve correction coefficients
%         optionsSelection(l).coeff = squeeze(coeff(k, min(find(strcmp(pixelValueSelectionMethods , optionsSelection(l).pixelValueSelectionMethod)),3), :))';
        [est(:,l), rmse(l,k), nmse(l,k)] = reflectanceEstimation( g, measured, ID(k), optionsSelection(l)); 
     end
    if (options.showImages)
        plots( 'estimationComparison', 2, [measured, est], sampleName, 'wavelength', wavelength, 'method',  options.pixelValueSelectionMethod, ...
            'saveOptions', options.saveOptions, 'lineNames', {'MS center \lambda', 'Measured', 'Est-Markovian', 'Est-Xcorr all data', 'Est-Xcorr macbeth', 'Est-Xcorr all sample','Est-Xcorr cancer sample', 'Est-Xcorr fixed cancer sample'});  
    end

%     plots('singlemeasurement', 3, wavelength, [], [], measured)
end

%% Export errors
rmse = rmse(:,rmse(1,:) ~= 0);
nmse = nmse(:,nmse(1,:) ~= 0);
errorData = [mean(rmse, 2), max(rmse,[], 2), min(rmse, [], 2), mean(nmse, 2), max(nmse,[], 2), min(nmse, [], 2)];
errors = struct('avgrmse', mean(rmse, 2), 'minrmse', min(rmse,[], 2), 'maxrmse',max(rmse, [], 2), 'stdrmse', std(rmse,[], 2), 'avgnmse', mean(nmse, 2), 'minnmse', min(nmse, [], 2), 'maxnmse', max(nmse,[], 2), 'stdnmse', std(nmse,[], 2), 'options', optionsSelection );

options.saveOptions.plotName = generateName(options, 'override', [], [], 'ComparisonRmseMatrixAvgMinMax');
plots('methodErrors', 3, [], 'RmseMatrixAvgMinMax', 'errors', errors, 'saveOptions', options.saveOptions)

options.saveOptions.plotName = generateName(options, 'override', [], [], 'ComparisonNmseMatrixAvgMinMax');
plots('methodErrors', 4, [], 'NmseMatrixAvgMinMax', 'errors', errors, 'saveOptions', options.saveOptions)

errorfilename = generateName(options, 'override', [], [], 'ComparisonErrorsMatrixAvgMinMax.mat');
save(errorfilename, 'errors');

%End of reflectance estimation matrix comparison

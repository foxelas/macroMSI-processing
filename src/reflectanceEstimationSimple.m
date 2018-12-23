%% Reflectance estimation at a selected small area / point    

%% main 
options.smoothingMatrixMethod = 'KCor same malignancy'; 
options.pixelValueSelectionMethod = 'extended';

rmse = zeros(1, msiN); 
nmse = zeros(1, msiN); 
estimatedSpectrumStruct = struct('Name', {}, 'Index', [], 'Spectrum', []);

% range = 146:187; %; %1:msiN; [66, 165, 175, 236, 297]
for k = 1:msiN
    % Give a name
    [options.saveOptions.plotName, sampleName] = generateName(options, 'plot+save',  data(ID(k).Representative), ID(k));    

    % Retrieve MSI data 
     g = MSIStruct(k).MSI; 

    % Retrieve measured spectrum
    measured = interp1(380:780, measuredSpectrumStruct(k).Spectrum, wavelength, 'nearest');   
    
   %% Estimation
%     options.coeff = squeeze(coeff(k, 3, :))';     
    [ est, rmse(1,k), nmse(1,k)] = reflectanceEstimation( g, measured, ID(k), options); 
    estimatedSpectrumStruct(k) = struct('Name', sampleName, 'Index', MSIStruct(k).Index, 'Spectrum', est);
    ID(k).Rmse = rmse(1,k);
    
    if (options.showImages)
        plots( 'estimationComparison', 2, [measured, est], sampleName, 'wavelength', wavelength, 'method', options.pixelValueSelectionMethod, ...
            'saveOptions', options.saveOptions, 'lineNames', {'f_c', 'Measured', 'Estimated'});  
    end
    
end

%% Export errors 
rmse = rmse(:,rmse(1,:) ~= 0);
nmse = nmse(:,nmse(1,:) ~= 0);
errors = struct('avgrmse', mean(rmse, 2), 'minrmse', min(rmse,[], 2), 'maxrmse',max(rmse, [], 2), 'stdrmse', std(rmse,[], 2), 'avgnmse', mean(nmse, 2), 'minnmse', min(nmse, [], 2), 'maxnmse', max(nmse,[], 2), 'stdnmse', std(nmse,[], 2), 'options', options );
out.EstimatedSpectrumStruct = estimatedSpectrumStruct;

% end of Reflectance estimation at a selected small area / point   
        
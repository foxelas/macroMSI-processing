%% Reflectance estimation of opposite cancer type 

%% print initial messages
disp('What happens if we use the cancerous correlation matrix for normal data and vice verca?')
disp('Were the doctors wrong into labeling skin samples during macropathology?')

%% main 
options.smoothingMatrixMethod = 'KCor all same malignancy'; 
options.pixelValueSelectionMethod = 'extended';

estimatedSpectrumStruct = struct('Name', {}, 'Index', [], 'Spectrum', []);
[~, ~, idx, ~] = subset('measured', name, 'unique');
count = 0;
falsePositive = 0;
falseNegative = 0;

% range = 146:187; %; %1:msiN; [66, 165, 175, 236, 297]
for j = 1:length(idx)
    k = idx(j);
    ID(k).IsNormal = ~ID(k).IsNormal;  

    % Retrieve MSI data 
     g = MSIStruct(k).MSI; 

    % Retrieve measured spectrum
    measured = interp1(380:780, measuredSpectrumStruct(k).Spectrum, wavelength, 'nearest');   
    
   %% Estimation
%     options.coeff = squeeze(coeff(k, 3, :))';     
    [ ~, rmse, nmse] = reflectanceEstimation( g, measured, ID(k), options);     
    
    ID(k).IsNormal = ~ID(k).IsNormal;

    if ID(k).Rmse - rmse > 0.0001 
        count = count + 1;
        if ID(k).IsNormal
            falsePositive = falsePositive + 1;
        else 
            falseNegative = falseNegative +1;
        end
    end 
    
end

%% Export errors 
disp('How many samples could have been labeled wrongly?')
fprintf('\nThe answer is: %d samples or %.2f%% out of the total measurements.\n', count, count / length(idx) * 100);
fprintf('False positives (labeled as normal but is cancer): %d(%.2f%%).\n', falsePositive, falsePositive / length(idx) * 100);
fprintf('False negatives (labeled as cancer but is normal): %d(%.2f%%).\n', falseNegative, falseNegative / length(idx) * 100);
% end of Reflectance estimation of opposite cancer type   
        
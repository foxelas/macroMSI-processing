%Start of reflectance estimation comparison
options.action = 'Refest';

w = warning('off', 'all');

methods = {'RGB-Simple', 'MSI-SpatioSpectral', 'MSI-Simple'}; 
methodsN = length(methods); 
estimatedSpectra = zeros(msiN, length(wavelength), methodsN);
nmses = zeros(msiN, methodsN);

for k = 1:msiN
    % Retrieve MSI data
    
    infile = fullfile(options.systemdir, 'infiles', strcat('poi_', num2str(k), '.mat'));
    load(infile, 'poiName', 'poiRAW', 'poiSegmentMask', ...
        'roiSeeds', 'measuredSpectrum', 'poiWhite');
        
    msi = poiRAW;
    mask = poiSegmentMask;
    measured = measuredSpectrum;
	
	[reconstructedArray, ~, nmseArray] = getMumtlipleReflectanceReconstructions( msi, rgb, mask, measured, ID(k), options, methods)
    estimatedSpectra(k,:,:) = reconstructedArray; 
	nmses(k,:) = nmseArray;
	
end
warning(w);

%% Export results
errorInfo = GetErrorInfoStruct(nmses);

for i = 1:methodsN
	nmseBars = getNmseForBarPlot(nmses(:, i), ID, methods{i}, options);
end

filename = mkdir_custom(fullfile(options.saveOptions.savedir, '5-ReflectanceEstimation', 'refest.mat'));
msiId = find(strcmp(methods, 'MSI-Simple'));
EstimatedSpectra =  reconstructedArray(:,:, msiId); 
rgbId = find(strcmp(methods, 'RGB-Simple'));
EstimatedRGBSpectra = reconstructedArray(:,:, rgbId); 
save(filename, 'EstimatedSpectra', 'EstimatedRGBSpectra','-append');
	
%End of reflectance estimation comparison


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function errorInfo = GetErrorInfoStruct(nmse, pixelValueSelectionMethods)
    errorInfo = struct('avgnmse', num2cell(mean(nmse, 2)), 'minnmse', num2cell(min(nmse, [], 2)), ...
                'maxnmse', num2cell(max(nmse, [], 2)), 'stdnmse', num2cell(std(nmse, [], 2)), 'stdeSEM',  num2cell(std(rmse, [], 2)./size(rmse,2)));
end


function nmseBars = getNmseForBarPlot(nmses, ID, method, options)
    i = 0;
    nmseBars = zeros(6, 1);
    for state = {'unfixed', 'fixed', 'cut'}
        for malignancy = 0:1
            i = i + 1;
			idx = [ID.IsBenign] ~= malignancy & strcmp({ID.Type};
            nmseBars(i) = mean(nmses(idx));
        end
        %fprintf('NRMSE for %s = %.4f\n', state{1}, nrmse);
    end
    nmseBars = reshape(nmseBars, [3,2])';
	
	options.saveOptions.plotName = fullfile(options.saveOptions.savedir, '6-ReflectanceEstimationPerformance', method);
    plotReconstructionPerformanceBars(nmseBars,{'benign', 'malignant'},'',1,options.saveOptions);
	
    fprintf('NRMSE overall = %.4f\n', mean(nmseBars(:)));
	
end



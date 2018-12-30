%% Reflectance estimation iterratively until miniminum estimation error point is reached

searchStep = 10;
width = 5;
height = 5;
options.smoothingMatrixMethod = 'KCor same malignancy';
options.pixelValueSelectionMethod = 'extended';

acceptedCount = 0;

for k = 1:msiN
    fprintf('Now optimization at %d\n', k);
    
    %% set options for min search
    sampleName = generateName([], 'image', data(ID(k).Representative), ID(k));
    
    measured = interp1(380:780, measuredSpectrumStruct(k).Spectrum, wavelength, 'nearest');
    %         options.coeff = squeeze(coeff(k, 3, :))';
    gOld = readMSI({data(ID(k).Data).File}, ID(k).Originx, ID(k).Originy, width, height);
    [estimatedOld, rmseOld, ~] = reflectanceEstimation(gOld, measured, ID(k), options);
    
    %finderr = @(xx, datafile, meas, id, sens, illum, opt) reflectanceEstimation( readMSI( datafile, fix(xx(1)), fix(xx(2)), 5, 5),  meas, id, sens, illum, opt);
    minfun = @(xx)finderror(xx, {data(ID(k).Data).File}, measured, ID(k), options);
    
    optionsSearch = optimset('MaxIter', 100);
    [xx, fval, ~, ~] = fminsearch(minfun, [ID(k).Originx, ID(k).Originy], optionsSearch);
    minLocation = fix(xx);
    
    %% To check before and after estimation
    gNew = readMSI({data(ID(k).Data).File}, minLocation(1), minLocation(2), width, height);
    [estimatedNew, rmseNew, ~] = reflectanceEstimation(gNew, measured, ID(k), options);
    
    if (options.showImages)
        plots('estimationComparison', 2, [measured, estimatedOld, estimatedNew], sampleName, 'wavelength', wavelength, 'method', options.pixelValueSelectionMethod, ...
            'lineNames', {'f_c', 'Measured', 'Old estimate', 'New Estimate'});
    end
    
    if (abs(ID(k).Originx-minLocation(1)) < 30 && abs(ID(k).Originy-minLocation(2)) < 30)
        ID(k).Originx = minLocation(1);
        ID(k).Originy = minLocation(2);
        ID(k).Rmse = rmseNew;
        acceptedCount = acceptedCount + 1;
        % disp('accepted')
    end
    
end
save([options.systemdir, 'IDmin.mat'], 'ID')
fprintf('Origin coordinates were adjusted for %d out of total %d samples', acceptedCount, msiN);

% %% To check new points on Image
ll = load('..\MATLAB\Data\saitamav2\ID.mat');
ID2 = ll.ID;
for i = 1:length(ID)
    idx = find([data.UniqueCount] == ID(i).UniqueCount, 1);
    showCroppedSection([], data(idx).File, ID(i).Originx, ID(i).Originy, [], [], true, 'y');
    showCroppedSection([], data(idx).File, ID2(i).Originx, ID2(i).Originy, [], [], false, 'm');
    pause(0.5)
end
clear ll

%end of Reflectance estimation iterratively until miniminum estimation error point is reached

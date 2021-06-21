function [tform] = getRegistrationTransform(msi1, msi2, method)

if nargin < 3
    method = 'surf';
end

if ndims(msi1) ~= ndims(msi2)
    disp('The fixed and moving image have different dimensions.');
    return;
end

if ndims(msi1) == 2
    fixed = msi1;
    moving = msi2;
elseif ndims(msi1) == 3
    fixed = squeeze(msi1(1, :, :));
    moving = squeeze(msi2(1, :, :));
elseif ndims(msi1) == 4
    fixed = squeeze(msi1(1, :, :, 1));
    moving = squeeze(msi2(1, :, :, 1));
end

savedir = getSetting('savedir');
childDir = getSetting('registration');
savedir = fullfile(savedir, childDir, method);

switch method
    case 'surf'
        ptsFixed = detectSURFFeatures(fixed);
        ptsMoving = detectSURFFeatures(moving);
        %             cornersFixed = detectFASTFeatures(fixed);
        %             cornersMoving= detectFASTFeatures(moving);
        [featuresFixed, validPtsFixed] = extractFeatures(fixed, ptsFixed);
        [featuresMoving, validPtsMoving] = extractFeatures(moving, ptsMoving);
        
        indexPairs = matchFeatures(featuresFixed, featuresMoving);
        matchedFixed = validPtsFixed(indexPairs(:, 1));
        matchedMoving = validPtsMoving(indexPairs(:, 2));
        
        setSetting('plotName', fullfile(savedir, strcat('fused1.png')));
        plots(1, @plotFused, fixed, moving, matchedFixed, matchedMoving, 'Putatively matched points (including outliers)');
        [tform, inlierMoving, inlierFixed] = estimateGeometricTransform( ...
            matchedMoving, matchedFixed, 'similarity');
        
        setSetting('plotName', fullfile(savedir, strcat('fused2.png')));
        plots(2, @plotFused, fixed, moving, inlierFixed, inlierMoving, 'Matching points (inliers only)');
        
    case 'regconfig'
        [optimizer, metric] = imregconfig('multimodal');
        movingRegisteredDefault = imregister(moving, fixed, 'affine', optimizer, metric);
        
        setSetting('plotName', fullfile(savedir, strcat('fused1.png')));
        plots(1, @plotMontage, fixed, movingRegisteredDefault, 'A: Default Registration');
        
        disp(optimizer)
        disp(metric)
        optimizer.InitialRadius = optimizer.InitialRadius / 3.5;
        movingRegisteredAdjustedInitialRadius = imregister(moving, fixed, 'affine', optimizer, metric);
        
        setSetting('plotName', fullfile(savedir, strcat('fused2.png')));
        plots(2, @plotMontage, fixed, movingRegisteredAdjustedInitialRadius, 'B: Adjusted InitialRadius');
        
        optimizer.MaximumIterations = 300;
        movingRegisteredAdjustedInitialRadius300 = imregister(moving, fixed, 'affine', optimizer, metric);
        setSetting('plotName', fullfile(savedir, strcat('fused3.png')));
        plots(3, @plotMontage, fixed, movingRegisteredAdjustedInitialRadius300, 'C: Adjusted InitialRadius, MaximumIterations = 300');
        
        tform = imregtform(moving, fixed, 'similarity', optimizer, metric);
        Rfixed = imref2d(size(fixed));
        movingRegisteredRigid = imwarp(moving, tform, 'OutputView', Rfixed);
        
        setSetting('plotName', fullfile(savedir, strcat('fused4.png')));
        plots(4, @plotMontage, fixed, movingRegisteredRigid, 'D: Registration Based on Similarity Transformation Model');
        
        movingRegisteredAffineWithIC = imregister(moving, fixed, 'affine', optimizer, metric, ...
            'InitialTransformation', tform);
        setSetting('plotName', fullfile(savedir, strcat('fused5.png')));
        plots(5, @plotMontage, fixed, movingRegisteredAffineWithIC, 'E: Registration from Affine Model Based on Similarity Initial Condition');
        
    otherwise
        disp('Unsupported registration method.');
end

Tinv = tform.invert.T;

ss = Tinv(2, 1);
sc = Tinv(1, 1);
scaleRecovered = sqrt(ss*ss+sc*sc)
thetaRecovered = atan2(ss, sc) * 180 / pi

outputView = imref2d(size(fixed));
recovered = imwarp(moving, tform, 'OutputView', outputView);

setSetting('plotName', fullfile(savedir, strcat('finalRegistered.png')));
plots(10, @plotMontage, fixed, recovered, 'Fixed vs Registered Moving');
end
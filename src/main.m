%% Set-up of options for running 
setup;

%% main 
tic;
fprintf('Running action %s...\n', options.action);    

switch (options.action)
    case 'ConfPlots'   
        %% Configuration plots about illumination, sensitivity etc
        %Checked
        options.saveOptions.plotName = strcat(options.saveOptions.savedir, 'illuminationPlot.jpg');
        plots( 'illumination', 1, [], '', 'wavelength', wavelength, 'illumination', illumination, 'saveOptions', options.saveOptions);
        options.saveOptions.plotName = strcat(options.saveOptions.savedir, 'sensitivityPlot.jpg');
        plots( 'sensitivity', 2, [], '', 'wavelength', wavelength, 'sensitivity', sensitivity, 'saveOptions', options.saveOptions);
        options.saveOptions.plotName = strcat(options.saveOptions.savedir, 'illuminationAndSensitivityPlot.jpg');
        plots( 'illuminationAndSensitivity', 3, [], '', 'wavelength', wavelength, 'illumination', illumination, 'sensitivity', sensitivity, 'saveOptions', options.saveOptions);
        % end of Configuration plots about illumination, sensitivity etc
    case 'CountData'
        %% Count and analyze the contents of the dataset
        [~, idx, ~] = unique(strcat({ID.Csvid}, {ID.T}));
        measuredSpectraCount = length(idx);
        normalFixedCount = 0;
        normalCutCount = 0;
        normalUnfixedCount = 0;
        fixedCount = 0;
        cutCount = 0;
        unfixedCount = 0;
        for i = 1:length(idx)
            if ~(ID(idx(i)).IsFixed)
                unfixedCount = unfixedCount + 1;
                if ID(idx(i)).IsNormal
                    normalUnfixedCount = normalUnfixedCount + 1;
                end
            elseif (ID(idx(i)).IsFixed && ~(ID(idx(i)).IsCut))
                fixedCount = fixedCount + 1;
                if ID(idx(i)).IsNormal
                    normalFixedCount = normalFixedCount + 1;
                end
            else
                cutCount = cutCount + 1;
                if ID(idx(i)).IsNormal
                    normalCutCount = normalCutCount + 1;
                end
            end           
        end
        normalCount = normalUnfixedCount+normalFixedCount+normalCutCount;
        cancerCount = length(idx) - normalCount;
        cancerUnfixedCount = unfixedCount - normalUnfixedCount;
        fprintf('Breakdown of the dataset:\n')
        fprintf('Total: %d, Unfixed: %d, Fixed: %d, Cut: %d\n', length(idx), unfixedCount, fixedCount, cutCount);
        fprintf('Normal: %d, Cancerous: %d\n\n', normalCount, cancerCount);
        fprintf('Among unfixed data:\n')
        fprintf('Normal: %d, Cancerous: %d\n\n', normalUnfixedCount, unfixedCount - normalUnfixedCount);
        fprintf('Among fixed data:\n')
        fprintf('Normal: %d, Cancerous: %d\n\n', normalFixedCount, fixedCount - normalFixedCount);
        fprintf('Among cut data:\n')
        fprintf('Normal: %d, Cancerous: %d\n\n', normalCutCount, cutCount - normalFixedCount);

%         Breakdown of the dataset:
%         Unfixed: 43, Fixed: 45, Cut: 39
%         Normal: 74, Cancerous: 53
% 
%         Among unfixed data:
%         Normal: 23, Cancerous: 20

        % end of  Count and analyze the contents of the dataset
    case 'FixMinimumErrorPoints'
        findMinimumRefErrorPoints;

    case 'ReflectanceEstimationSystemComparison'  
        reflectanceEstimationSystemComparison;
        
    case 'ReflectanceEstimationMatrixComparison'
        reflectanceEstimationMatrixComparison;
        
    case 'ReflectanceEstimationMatrixSystemComparison'
        reflectanceEstimationMatrixSystemComparison;
        
    case 'ReflectanceEstimationNoiseComparison'
        reflectanceEstimationNoiseComparison;
        
    case 'ReflectanceEstimationSimple'
        reflectanceEstimationSimple;
    
    case 'ReflectanceEstimationOpposite'
        reflectanceEstimationOpposite;
        
    case 'CreateSRGB'  % from the MSI reflectances 
        %% CREATE SRGB
        
        for k = [ 250, 270 ]
            [~, sampleName] = generateName(options, 'plot+save',  data(ID(k).Representative), ID(k));    
             g = readMSI({data(ID(k).Data).File}, [],[],[],[], [0, 450, 465, 505, 525, 575, 605, 630 ], true ); 
            sRGB = createSRGB(g, 'original', sampleName, ID(k), options);             
        end  
        
        % end of create srgb
        
    case 'PCA' 
        %% PCA for spectrum data
        input = 'estimated';
        [Gun, lineNamesun] = subset(input, name, 'unique');        
        % G 401 variables x 39 observations 
        %[coeff, score, latent, ~, explained]= pca(Gun, 'Centered', true, 'Economy', false);   % has inversion error 
        [~, score, latent, explained] = PCA(Gun);
        options.saveOptions.plotName = generateName(options, 'override', [], [], ['PCA of ', input, ' spectra by type']);
        plots('pca', 1, score, 'PCA Fix', 'lineNames', lineNamesun, 'latent', latent(1:10), 'explained', explained(1:10), 'saveOptions', options.saveOptions)
        options.saveOptions.plotName = generateName(options, 'override', [], [], ['PCA of ', input, ' spectra by sample']);
        plots('pca', 2, score, 'PCA Sample', 'lineNames', lineNamesun, 'latent', latent(1:10), 'explained', explained(1:10), 'saveOptions', options.saveOptions)
        
        %unfixed data only 
        [Gfx, lineNamesfx] = subset(input, name, 'unfixed');          
        %[coeff, score, latent, ~, explained]= pca(Gfx,'Centered', true);   
        [~, score, latent, explained] = PCA(Gfx);
        options.saveOptions.plotName = generateName(options, 'override', [], [], ['PCA of ', input,  'spectra (Only Unfixed case)']);
        plots('pca', 3, score, 'PCA Sample', 'lineNames', lineNamesfx, 'latent', latent(1:10), 'explained', explained(1:10), 'saveOptions', options.saveOptions)
        
        input = 'measured';
        [Gun, lineNamesun] = subset(input, name, 'unique');          
        %[coeff, score, latent, ~, explained]= pca(Gun,'Centered', true);   
        [~, score, latent, explained] = PCA(Gun);
        options.saveOptions.plotName = generateName(options, 'override', [], [], ['PCA of ', input, ' spectra by type']);
        plots('pca', 4, score, 'PCA Fix', 'lineNames', lineNamesun, 'latent', latent(1:10), 'explained', explained(1:10), 'saveOptions', options.saveOptions)
        options.saveOptions.plotName = generateName(options, 'override', [], [], ['PCA of ', input, ' spectra by sample']);
        plots('pca', 5, score, 'PCA Sample', 'lineNames', lineNamesun, 'latent', latent(1:10), 'explained', explained(1:10), 'saveOptions', options.saveOptions)
    
        %unfixed data only 
        [Gfx, lineNamesfx] = subset(input, name, 'unfixed');          
        [coeff, score, latent, explained] = PCA(Gfx);
        options.saveOptions.plotName = generateName(options, 'override', [], [], ['PCA of ', input, ' spectra (Only Unfixed case)']);
        plots('pca', 6, score, 'PCA Sample', 'lineNames', lineNamesfx, 'latent', latent(1:10), 'explained', explained(1:10), 'saveOptions', options.saveOptions)
       
        % End PCA for spectrum data      
        
    case 'LDA'
         %% LDA for spectrum data 
        input = 'estimated';
        ldaMethod = 'p'; % 'pdf', 'kardi', ''
        [Gun, lineNamesun, ~, labelsun] = subset(input, name, 'unique');          
        [~, scores] = LDA(Gun, double(labelsun), [], ldaMethod);
        %L = [ones(length(un),1) Gun] * W';
        %P = exp(L) ./ repmat(sum(exp(L),2),[1 2]);
        options.saveOptions.plotName = generateName(options, 'override', [], [], ['LDA of ', input, ' spectra by type', ldaMethod]);
        plots('lda', 1, scores, 'LDA Fix', 'lineNames', lineNamesun, 'saveOptions', options.saveOptions)
        options.saveOptions.plotName = generateName(options, 'override', [], [], ['LDA of ', input, ' spectra by sample', ldaMethod]);
        plots('lda', 2, scores, 'LDA Sample', 'lineNames', lineNamesun, 'saveOptions', options.saveOptions)

        %unfixed data only 
        [Gfx, lineNamesfx, ~, labelsfx] = subset('estimated', name, 'unfixed');          
        [~, scores] = LDA(Gfx, double(labelsfx), [], ldaMethod);
        options.saveOptions.plotName = generateName(options, 'override', [], [], ['LDA of ', input, ' spectra by sample (only unfixed)', ldaMethod]);
        plots('lda', 3, scores, 'LDA Sample', 'lineNames', lineNamesfx, 'saveOptions', options.saveOptions)

        input = 'measured';
        [Gun, lineNamesun, ~, labelsun] = subset(input, name, 'unique');          
        [~, scores] = LDA(Gun, double(labelsun), [], ldaMethod);
        %L = [ones(length(un),1) Gun] * W';
        %P = exp(L) ./ repmat(sum(exp(L),2),[1 2]);
        options.saveOptions.plotName = generateName(options, 'override', [], [], ['LDA of ', input, ' spectra by type', ldaMethod]);
        plots('lda', 4, scores, 'LDA Fix', 'lineNames', lineNamesun, 'saveOptions', options.saveOptions)
        options.saveOptions.plotName = generateName(options, 'override', [], [], ['LDA of ', input, ' spectra by sample', ldaMethod]);
        plots('lda', 5, scores, 'LDA Sample', 'lineNames', lineNamesun, 'saveOptions', options.saveOptions)

        %unfixed data only 
        [Gfx, lineNamesfx, ~, labelsfx] = subset('estimated', 'unfixed', options.saveOptions.savedir);          
        [~, scores] = LDA(Gfx, double(labelsfx), [], ldaMethod);
        options.saveOptions.plotName = generateName(options, 'override', [], [], ['LDA of ', input, ' spectra by sample (only unfixed)', ldaMethod]);
        plots('lda', 6, scores, 'LDA Sample', 'lineNames', lineNamesfx, 'saveOptions', options.saveOptions)

         % End LDA for spectrum data 
    case 'PCA-LDA'
        %% Do first PCA and then LDA
        input = 'estimated';
        [Gfx, lineNamesfx, fx, labelsfx] = subset(input, name, 'unfixed');          

        %Lower dimensions to 10 
        [~, scores, latent, ~, explained]= pca(Gfx,'Centered', true);   
        options.saveOptions.plotName = generateName(options, 'override', [], [], ['PCA of ', input, ' spectra (Only Unfixed case)']);
        plots('pca', 1, scores, 'PCA Sample', 'lineNames', lineNamesfx, 'latent', latent(1:10), 'explained', explained(1:10), 'saveOptions', options.saveOptions)
        
        %Class separation 
        ldaMethod = 'pdf';
        [~, scores] = LDA(scores(:, 1:10), double(labelsfx), [], ldaMethod);
        options.saveOptions.plotName = generateName(options, 'override', [], [], ['LDA after PCA of ', input, ' spectra by sample (only unfixed)', ldaMethod]);
        plots('lda', 2, scores, 'LDA Sample', 'lineNames', lineNamesfx, 'saveOptions', options.saveOptions) 

        % End first PCA and then LDA
        
    case 'Nearest'
        version = 'estimated';
        
        classificationError = struct('Accuracy', [], 'TypeI', [], 'TypeII', [], 'Input', {}, 'Projection',{}, 'Validation', {}, 'VoteRule', {}, 'Neighbours', [], 'Distance', {});
        validation = {'LeaveMOut', 'Kfold'};
        validationShort = {'1O', 'Kf'};
        nnrRules = {'majority', 'weighted majority', 'complex vote'};
        nnrRulesShort = {'M', 'W', 'C'};
        distances = {'correlation', 'chebychev', 'euclidean'};
        distancesShort = {'Cr', 'Ch', 'Eu'};
        groups = {'unique', 'fixed', 'unfixed'};
        projections = {'PCA', 'LDA', 'PCALDA', 'Spectrum'};
        options.saveOptions.saveInHQ = true;
        for j = 1:length(validation)
            classificationError = struct('Accuracy', [], 'TypeI', [], 'TypeII', [], 'Input', {}, 'Projection',{}, 'Validation', {}, 'VoteRule', {}, 'Neighbours', [], 'Distance', {});
            m = 0;
            for g = 1:length(groups)
                for p = 1:length(projections)
                    [Gun, labels] = classifierInput(version, groups{g}, projections{p}, name);
                    for i = 1:length(nnrRules)
                        for k = [1 3 5]
                            for d = 1:length(distances)
                                [a, b, c] = knn(Gun, labels, k, distances{d}, nnrRules{i}, validation{j}); % labels = 1 for cancer 
                                m = m + 1;
                                classificationError(m) = struct('Accuracy', a, 'TypeI', b, 'TypeII', c, 'Input', groups{g}, 'Projection', projections{p}, ...
                                    'Validation', validation{j}, 'VoteRule', nnrRules{i}, 'Neighbours', k, 'Distance', distances{d});
                            end
                        end
                    end
                end               
            end
            classificationErrorFields = fieldnames(classificationError);
            classificationErrorCell = struct2cell(classificationError);
            sz = size(classificationErrorCell);
            % Convert to a matrix
            classificationErrorCell = reshape(classificationErrorCell, sz(1), []); 
            % Make each field a column
            classificationErrorCell = classificationErrorCell';
            % Sort by first field "name"
            classificationErrorCell = sortrows(classificationErrorCell, -1);
            classificationErrorCell = reshape(classificationErrorCell', sz);
            classificationError = cell2struct(classificationErrorCell, classificationErrorFields, 1);
            options.saveOptions.plotName = generateName(options, 'override', [], [], ['Classification error of ', version, ' spectra with ', validation{j}]);
            plots('classificationErrors', 2, [], '', 'errors', classificationError, 'saveOptions', options.saveOptions)
        end

                
    case {'MeasuredOverlap', 'EstimatedOverlap'}
        %% Overlap of cancer / fixed
        
        if containts(options.action, 'Measured')
            [~, xx] = unique(strcat({ID.Csvid}, {ID.T}), 'last');
            id = ID(xx);
            curves = zeros( length(380:780), length(id));
        else 
            e = matfile(fullfile(options.saveOptions.saveDir, 'ReflectanceEstimationSimple', 'out.mat')); 
            id = ID;
            curves = e.EstimatedSpectrumStruct;
        end
        
        %cancer
        lineNamesCancer = cell( length(id), 1);        
        lineNamesFixed = cell( length(id), 1);        
        lineNamesSample = cell( length(id), 1);        

        for k = 1:length(id)
            curves(:,k) = measuredSpectrumStruct(k).Spectrum;
            if id(k).IsNormal
                lineNamesCancer{k} = 'normal';
            else
                lineNamesCancer{k} = 'cancerous';
            end
            if id(k).IsFixed
                lineNamesFixed{k} = 'fixed';
            else
                lineNamesFixed{k} = 'unfixed';
            end
            if id(k).IsNormal
                lineNamesSample{k} = [id(k).Sample, '_normal'];
            else
                lineNamesSample{k} = [id(k).Sample, '_cancerous'];
            end
        end
        options.saveOptions.plotName = generateName(options, 'override', [], [], [ options.action, 'Cancer']);
        plots( 'overlapSpectrum', 2, 380:780, curves, 'Cancer', 'wavelength', wavelength, 'lineNames', lineNamesCancer, 'saveOptions', options.saveOptions);
        options.saveOptions.plotName = generateName(options, 'override', [], [], [ options.action, 'Fixed']);
        plots( 'overlapSpectrum', 3, 380:780, curves, 'Fixed', 'wavelength', wavelength, 'lineNames', lineNamesFixed, 'saveOptions', options.saveOptions);
        options.saveOptions.plotName = generateName(options, 'override', [], [], [ options.action, 'Sample']);
        plots( 'overlapSpectrum', 1, 380:780, curves, 'Sample', 'wavelength', wavelength, 'lineNames', lineNamesSample, 'saveOptions', options.saveOptions);
        
        % end overlap of cancer / fixed
        
    otherwise
        disp('Nothing to do here. Aborting...');
end
fprintf('Finished running action %s.\n', options.action);    
elapsed = toc;

outputLog = logInfo(options, elapsed);

% actionING ends



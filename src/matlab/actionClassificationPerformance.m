currdir = 'saitama_v8_min_region_bright';
filename = fullfile('..', '..', '..', 'output', currdir, 'msiclas.csv');
%filename = fullfile('..', '..', '..', 'output', currdir, 'rgbclas.csv');
%filename = fullfile('..', '..', '..', 'output', currdir, 'classification_logs', 'RGB_validation_2scales2019-07-25 23_25.csv');

fixing = { 'unfixed','fixed', 'mixed'};
if contains(filename, 'only_spect')
    savedir = fullfile('..', '..', '..', 'output', currdir, '10-ClassifierPerformanceOnlySpect' );
    lbps = {'spect'};
    classifiers = {'SVM', 'KNN'}; % , 'LDA' , 'QDA'
else
    savedir = fullfile('..', '..', '..', 'output', currdir, '9-ClassifierPerformance' );
    classifiers = {'SVM', 'KNN', 'Random Forest'};
end

saveOptions.savedir = savedir;
saveOptions.saveImages = true;
saveOptions.plotName = '';
saveOptions.saveInHQ = false;
[tab, dimred, class_s] = ReadClassificationData(filename);
summary(tab)

if contains(lower(filename), 'rgb')
    currentCase = 'RGB-reconstructed Spectra';
    saveOptions.savedir = strcat(saveOptions.savedir, '_rgb');
    %lbps = {'spect', 'spect+SumLBP'};
    lbps = {'spect', 'spect+MMLBP'};
    lbpsShow =  {'spect', 'spect+LBP'};
elseif contains(filename, 'mes')
    currentCase = 'Measured Spectra';
else 
    currentCase = 'MSI-reconstructed Spectra';
    lbps = {'spect','spect+CatLBP', 'spect+MMLBP', 'spect+SumLBP'};
    lbpsShow = lbps;
end

if contains(filename, 'only_spect')
    cond1 = classifiers;
    cond2 = {'fixed', 'unfixed', 'mixed'};
    [barAccInfo1, maxIds1] = performanceComparisonArray(tab, class_s, cond1, cond2, 'Specificity', 'BalancedAccuracy');
    barAccInfo2 = reshape(tab(reshape(maxIds1, [1,6]), 'BalancedAccuracy').BalancedAccuracy, [2,3]);
    plots('classificationPerformanceBars', 2, [], {strcat(['Performance for ', currentCase]), 'Accuracy', 'Specificity'}, 'LineNames', ...
        classifiers, 'Performance', {barAccInfo2, barAccInfo1}, 'SaveOptions', saveOptions);
else
   
    disp('Classifier comparison')
    saveOptions.plotName = fullfile(saveOptions.savedir, 'unfixed_classifier_auc');
    [barAccInfo1, idxs] = performanceComparisonArray(tab, class_s, lbps,  classifiers , 'AUC', 'Accuracy', 'unfixed');
    plotClassificationPerformanceBars(barAccInfo1, lbpsShow, {'SVM', 'KNN', 'RF'}, 'ROC AUC', [50, 100], 1, saveOptions);
    saveOptions.plotName = fullfile(saveOptions.savedir, 'unfixed_classifier_dor');
    barDORInfo1 = getRespectiveValues(barAccInfo1, idxs, tab.DOR);
    plotClassificationPerformanceBars(barDORInfo1 ./ 100, lbpsShow, {'SVM', 'KNN', 'RF'}, 'DOR', [0, 30], 1, saveOptions);
    
    saveOptions.plotName = fullfile(saveOptions.savedir, 'mixed_classifier_auc');
    [barAccInfo1, idxs] = performanceComparisonArray(tab, class_s,   lbps,  classifiers , 'AUC', 'Accuracy', 'mixed');
    plotClassificationPerformanceBars(barAccInfo1, lbpsShow, {'SVM', 'KNN', 'RF'}, 'ROC AUC', [50, 100], 2, saveOptions);
    
    saveOptions.plotName = fullfile(saveOptions.savedir, 'mixed_classifier_dor');
    barDORInfo1 = getRespectiveValues(barAccInfo1, idxs, tab.DOR);
    plotClassificationPerformanceBars(barDORInfo1 ./ 100, lbpsShow, {'SVM', 'KNN', 'RF'}, 'DOR', [0, 30], 2, saveOptions);
    
    disp('Fixing comparison')
    [barAccInfo2, idx1] = performanceComparisonArray(tab, class_s, lbps, fixing, 'AUC', 'Accuracy', 'KNN');
    saveOptions.plotName = fullfile(saveOptions.savedir, 'knn_tissue_auc');
    plotClassificationPerformanceBars(barAccInfo2, lbpsShow, {'Unfixed', 'Fixed', 'Mixed'}, 'ROC AUC', [50, 100], 3, saveOptions);
    [barAccInfo2, idx2] = performanceComparisonArray(tab, class_s, lbps, fixing, 'AUC', 'Accuracy', 'SVM');
    saveOptions.plotName = fullfile(saveOptions.savedir, 'svm_tissue_auc');
    plotClassificationPerformanceBars(barAccInfo2, lbpsShow, {'Unfixed', 'Fixed', 'Mixed'},  'ROC AUC', [50, 100], 4, saveOptions);
    [barAccInfo2, idx3] = performanceComparisonArray(tab, class_s, lbps, fixing, 'AUC', 'Accuracy', 'Random Forest');
    saveOptions.plotName = fullfile(saveOptions.savedir, 'rf_tissue_auc');
    plotClassificationPerformanceBars(barAccInfo2, lbpsShow, {'Unfixed', 'Fixed', 'Mixed'},  'ROC AUC', [50, 100], 5, saveOptions);

    barDORInfo1 = getRespectiveValues(barAccInfo1, idx1, tab.DOR);
    saveOptions.plotName = fullfile(saveOptions.savedir, 'knn_tissue_dor');
    plotClassificationPerformanceBars(barDORInfo1 ./ 100, lbps, {'Unfixed', 'Fixed', 'Mixed'}, 'DOR', [0, 30], 3, saveOptions);
    barDORInfo1 = getRespectiveValues(barAccInfo1, idx2, tab.DOR);
    saveOptions.plotName = fullfile(saveOptions.savedir, 'svm_tissue_dor');
    plotClassificationPerformanceBars(barDORInfo1 ./ 100, lbps, {'Unfixed', 'Fixed', 'Mixed'}, 'DOR', [0, 30], 4, saveOptions);
    barDORInfo1 = getRespectiveValues(barAccInfo1, idx3, tab.DOR);
    saveOptions.plotName = fullfile(saveOptions.savedir, 'rf_tissue_dor');
    plotClassificationPerformanceBars(barDORInfo1 ./ 100, lbps, {'Unfixed', 'Fixed', 'Mixed'}, 'DOR', [0, 30], 5, saveOptions);
    
    disp('Dimred comparison')
    %lbps = {'spect|fixed','spect+CatLBP|fixed', 'spect+MMLBP|fixed',
    %'spect+SumLBP|fixed'};
    cs = {'PCA10+None', 'PCA10+PCA', 'PCA10+ICA'};
    % cs = {'None+', 'PCA+', 'ICA+'};
    [barAccInfo2, ~] = performanceComparisonArray(tab, class_s, lbps, cs, 'AUC', 'Accuracy');
end 


%% Function for finding perfomance parameter for all conditions
function [barAccInfo, maxIds] = performanceComparisonArray(tab, class_s, cond1, cond2, comparison_param, additional_param, additional_cond)

    if (nargin < 7)
        additional_cond = [];
    end
    
    barAccInfo = zeros(numel(cond1),numel(cond2));
    maxIds = zeros(numel(cond1),numel(cond2));
    
    for ii = 1:numel(cond1)
        for jj = 1:numel(cond2)
            [barAccInfo(ii, jj), maxIds(ii, jj) ]= findMaxContaining(tab, comparison_param, additional_param, class_s, cond1{ii}, cond2{jj}, additional_cond);
        end
    end
    
end

function vals = getRespectiveValues(barAccInfo1, idxs, column)
    dorVals = zeros(size(barAccInfo1(:)));
    pos = find(idxs > 0);
    idxs2 = idxs(idxs > 0);
    dorVals(pos) = column(idxs2);
    vals = reshape( dorVals , size(barAccInfo1));
end

%% Function for Finding maximum AUC in performance data
function [maxComp, maxId] = findMaxContaining( tab, comparison_param, additional_param, class_s, condition1, condition2, condition3)

    if strcmp(condition1, 'spect')
        condition1 = 'spect|';
    end
    if contains(condition2, 'xed')
        condition2 = strcat('|', condition2);
    end
    
    if exist('condition3', 'var') && ~isempty(condition3)
        condIds = find(contains(class_s, condition1) & contains(class_s, condition2) & contains(class_s, condition3) ); %  & tab.AUC > 0.7 & tab.Accuracy > 0.6 & tab.AUC > 0.6);  
    else
        condIds = find(contains(class_s, condition1) & contains(class_s, condition2) & tab.AUC > 0.7); %  & tab.Accuracy > 0.6  & (~contains(condition2, 'SVM') | contains(classifier, 'sigmoid')) 
    end
    comparison_column = tab.(comparison_param);
    comparison_column = comparison_column(condIds);
    [maxComp, maxId] = max(comparison_column);
    additional_column = tab.(additional_param);
    additional_column = additional_column(condIds);
    maxAdd = additional_column(maxId);
    if isempty(condIds) && isempty(maxComp)
        maxAdd = 0;
        maxComp = 0;
        maxId = 0;
        fprintf(' : %s %.4f and %s %.4f \n', comparison_param, maxComp, additional_param, maxAdd);
    else
        maxId = condIds(maxId);
        fprintf('%s, %s: %s %.4f and %s %.4f \n', tab.Classifier{maxId},class_s{maxId}, comparison_param, maxComp, additional_param, maxAdd);
    end

end

%% Function for Reading performance data
function [tab, dimred, class_s] = ReadClassificationData(filename)
    tab = readtable(filename, 'Delimiter', ',', 'ReadVariableNames', true, 'Format','%s%s%s%s%s%s%s%f%f%f%f%f%f%f%f');
    %tab.Properties.VariableNames(1) = {'Classifier'};

    for i  = 1:length(tab.Classifier)
        if strcmp(tab.Input(i), 'unique_fixed')
            tab.Input{i} = 'fixed';
        elseif strcmp(tab.Input(i), 'unique_unfixed')
            tab.Input{i} = 'unfixed';
        else 
            tab.Input{i} = 'mixed';
        end
        tab.Feature{i} = strrep(tab.Feature{i}, 'clbp', 'CatLBP');
        tab.Feature{i} = strrep(tab.Feature{i}, 'slbp', 'SumLBP');
        tab.Feature{i} = strrep(tab.Feature{i}, 'mlbp', 'MMLBP');
        tab.DOR(i) = tab.Sensitivity(i) * tab.Specificity(i) / ((1 - tab.Sensitivity(i)) * (1 - tab.Specificity(i)));
    end

    dimred = cellfun(@(w,x,y,z) strjoin({strcat(w,num2str(x)),strcat(y, num2str(z))}, '+'),...
        tab.Dimred1, tab.NComp1, tab.Dimred2, tab.NComp2, 'uni', 0);

    class_s = cellfun(@(w,x,y,z) strrep(strjoin({w,x,y,z}, '|'), '_', ' ' ),...
        tab.Classifier, tab.Feature, tab.Input, dimred, 'uni', 0);
end
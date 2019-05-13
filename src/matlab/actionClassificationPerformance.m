%filename = fullfile('macroMSI-processing', 'logs', '2019-04-15 03_32.log');
%filename = fullfile('macroMSI-processing', 'logs', 'Classification_v2.log');
% filename = fullfile('macroMSI-processing', 'logs', 'Classification_v2_rgb.log');
filename = fullfile('..', '..', '..', 'output', 'saitama_v7_min_region_e', 'classification_only_spect.csv');

fixing = { 'unfixed','fixed', 'mixed'};
if contains(filename, 'only_spect')
    savedir = fullfile('..', '..', '..', 'output', 'saitama_v7_min_region_e', 'ClassifierPerformanceOnlySpect' );
    lbps = {'spect'};
    classifiers = {'SVM', 'KNN'}; % , 'LDA' , 'QDA'
else
    savedir = fullfile('..', '..', '..', 'output', 'saitama_v7_min_region_e', 'ClassifierPerformance' );
    lbps = {'spect','spect+CatLBP', 'spect+MMLBP', 'spect+SumLBP'};
    classifiers = {'SVM', 'KNN', 'RF'};
end
saveOptions.savedir = savedir;
saveOptions.saveImages = true;
saveOptions.plotName = '';
saveOptions.saveInHQ = false;
[tab, dimred, class_s] = ReadClassificationData(filename);
summary(tab)

if contains(filename, 'rgb')
    currentCase = 'RGB-reconstructed Spectra';
elseif contains(filename, 'mes')
    currentCase = 'Measured Spectra';
else 
    currentCase = 'MSI-reconstructed Spectra';
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
    [barAccInfo1, ~] = performanceComparisonArray(tab, class_s, lbps, classifiers, 'AUC', 'Accuracy');
    disp('Fixing comparison')
    [barAccInfo2, ~] = performanceComparisonArray(tab, class_s, lbps, fixing, 'AUC', 'Accuracy');
    plots('classificationPerformanceBars', 12, [], {strcat('Performance for ', currentCase), 'classifier', 'fixing'}, 'LineNames', lbps,...
        'Performance', {barAccInfo1, barAccInfo2}, 'SaveOptions', saveOptions);

    disp('Dimred comparison')
    %lbps = {'spect|fixed','spect+CatLBP|fixed', 'spect+MMLBP|fixed',
    %'spect+SumLBP|fixed'};
    cs = {'PCA10+None', 'PCA10+PCA', 'PCA10+ICA'};
    % cs = {'None+', 'PCA+', 'ICA+'};
    [barAccInfo2, ~] = performanceComparisonArray(tab, class_s, lbps, cs, 'AUC', 'Accuracy');
end 


%% Function for finding perfomance parameter for all conditions
function [barAccInfo, maxIds] = performanceComparisonArray(tab, class_s, cond1, cond2, comparison_param, additional_param)
    barAccInfo = zeros(numel(cond1),numel(cond2));
    maxIds = zeros(numel(cond1),numel(cond2));
    
    for ii = 1:numel(cond1)
        for jj = 1:numel(cond2)
            [barAccInfo(ii, jj), maxIds(ii, jj) ]= findMaxContaining(tab, comparison_param, additional_param, class_s, cond1{ii}, cond2{jj}, []);
        end
    end
    
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
        condIds = find(contains(class_s, condition1) & contains(class_s, condition2) & contains(class_s, condition3) & tab.Accuracy > 0.6 & tab.AUC > 0.6);  
    else
        condIds = find(contains(class_s, condition1) & contains(class_s, condition2)); %  & tab.Accuracy > 0.6 & tab.AUC > 0.6 & (~contains(condition2, 'SVM') | contains(classifier, 'sigmoid')) 
    end
    comparison_column = tab.(comparison_param);
    comparison_column = comparison_column(condIds);
    [maxComp, maxId] = max(comparison_column);
    additional_column = tab.(additional_param);
    additional_column = additional_column(condIds);
    maxAdd = additional_column(maxId);
    if isempty(maxComp)
        maxAdd = 0;
        maxComp = 0;
    end
    maxId = condIds(maxId);
    fprintf('%s, %s: %s %.4f and %s %.4f \n', tab.Classifier{maxId},class_s{maxId}, comparison_param, maxComp, additional_param, maxAdd);
end

%% Function for Reading performance data
function [tab, dimred, class_s] = ReadClassificationData(filename)
    tab = readtable(filename);
    tab.Properties.VariableNames(1) = {'Classifier'};

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
    end

    dimred = cellfun(@(w,x,y,z) strjoin({strcat(w,num2str(x)),strcat(y, num2str(z))}, '+'),...
        tab.Dimred1, tab.NComp1, tab.Dimred2, tab.NComp2, 'uni', 0);

    class_s = cellfun(@(w,x,y,z) strrep(strjoin({w,x,y,z}, '|'), '_', ' ' ),...
        tab.Classifier, tab.Feature, tab.Input, dimred, 'uni', 0);
end
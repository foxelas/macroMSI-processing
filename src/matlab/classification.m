savedir = fullfile('..', '..', '..', 'output', 'saitama_v7_min_region_e', 'ClassifierPerformance' );

filename = fullfile('macroMSI-processing', 'logs', '2019-04-15 03_32.log');
res = delimread(filename, ',', 'raw');
res = res.raw;
%filename = fullfile('..', '..', '..', 'output', 'saitama_v7_min_region_e', 'classification.xlsx');
%[num,txt,raw] = xlsread(filename, 'ACC');
classifier = res(:,1);
inputSet = res(:,2);
featureSet = res(:,3);
dim1 = res(:,4);
dim1num = res(:,5);
dim2 = res(:,6);
dim2num = res(:,7);
acc = str2double(res(:,8));
auc = str2double(res(:,9));

acceptedId = [];
for i  = 1:length(classifier)
    if strcmp(inputSet(i), 'unique_fixed')
        inputSet{i} = 'fixed';
    elseif strcmp(inputSet(i), 'unique_unfixed')
        inputSet{i} = 'unfixed';
    else 
        inputSet{i} = 'mixed';
    end
    featureSet{i} = strrep(featureSet{i}, 'clbp', 'CatLBP');
    featureSet{i} = strrep(featureSet{i}, 'slbp', 'SumLBP');
    featureSet{i} = strrep(featureSet{i}, 'mlbp', 'MMLBP');

    if ~(acc(i) < 0.75 || auc(i) < 0.75)
        acceptedId = [acceptedId, i];
    end
end

dimred = cellfun(@(w,x,y,z) strjoin({strcat(w),strcat(y)}, '+'),...
    dim1, dim1num, dim2, dim2num, 'uni', 0);

% figure(1);
% histogram(categorical(classifier(acceptedId)), 'Normalization', 'probability')
% title('Acceptable classifiers by Accuracy', 'FontSize', 15);
% xlabel('Classifiers')
% ylabel('Probability')
% 
% figure(2);
% histogram(categorical(featureSet(acceptedId)),'Normalization', 'probability')
% title('Acceptable classifiers by Features', 'FontSize', 15);
% xlabel('Feature Set')
% ylabel('Probability')
% 
% figure(3);
% histogram(categorical(inputSet(acceptedId)),'Normalization', 'probability')
% title('Acceptable classifiers by Input', 'FontSize', 15);
% xlabel('Feature Set')
% ylabel('Probability')
% 
% figure(4);
% histogram(categorical(dimred(acceptedId)), 'Normalization','probability')
% title('Acceptable classifiers by Dimension Reduction', 'FontSize', 15);
% xlabel('Dimension Reduction')
% ylabel('Probability')

markers = cell(length(classifier),1);
colors = cell(length(classifier),1);
for i  = 1:length(classifier)
    if contains(classifier(i), 'KNN')
        classifier{i} = 'KNN';
        markers{i} = '^';
        colors{i} = 'b';
    elseif contains(classifier(i), 'SVM')
        classifier{i} = 'SVM';
        markers{i} = 'x';
        colors{i} = 'm';
    else 
        classifier{i} = 'RF';
        markers{i} = 'o';
        colors{i} = 'g';
    end
end  
[acc_s, id1] = sort(acc, 'descend');
class_s = cellfun(@(w,x,y) strrep(strjoin({w,x,y}, '|'), '_', ' ' ),...
    classifier(id1), featureSet(id1), inputSet(id1), 'uni', 0);

auc_s = auc(id1);
figure(5);
yyaxis left
scatter([1:10],  acc_s(1:10));
ylim([0.75, 1])
ylabel('Accuracy')

yyaxis right
scatter([1:10], auc_s(1:10));
ylim([0.65, 1])
ylabel('Area Under Curve')
xticklabels(class_s);

title('Classifier Performance')
xtickangle(45);

plotClassificationPerformance(6, auc(acceptedId), acc(acceptedId), classifier(acceptedId), fullfile(savedir, 'Classifiers'));
plotClassificationPerformance(7, auc(acceptedId), acc(acceptedId), featureSet(acceptedId), fullfile(savedir, 'FeatureSet'));
plotClassificationPerformance(8, auc(acceptedId), acc(acceptedId), dimred(acceptedId), fullfile(savedir, 'DimRed'));
plotClassificationPerformance(9, auc(acceptedId), acc(acceptedId), inputSet(acceptedId), fullfile(savedir, 'InputSet'));
plotClassificationPerformance(10, auc(acceptedId), acc(acceptedId),  {classifier(acceptedId), featureSet(acceptedId) }, fullfile(savedir, 'Classifier+Features'));


function [] = plotClassificationPerformance(fig, aucVal, accVal, grouping, plotName)
figure(fig);
if (size(grouping,2) > 1) 
    gscatter(aucVal, accVal, grouping, 'rrrrcccc','.*^x+', [10, 10]);
else 
    gscatter(aucVal, accVal, grouping, 'rbgck', '.x', 10);
end
xlabel('AUC');
ylabel('Accuracy');
title('Perfomance of Top Classifiers', 'FontSize', 15)
legend('Location','northeastoutside')
ylim([0.7, 0.9])
xlim([0.7, 0.9])

if (~isempty(plotName))
    filename = strrep(plotName, '.mat', '');
    
    [filepath,name,~] = fileparts(filename);
    filepathBW = fullfile(filepath, 'bw');
    mkdir_custom(filepath);
    mkdir_custom(filepathBW);
   
    filename = fullfile(filepath, strcat(name, '.jpg'));
    export_fig(filename , '-jpg');
    filename = fullfile(filepath, strcat(name, '.eps'));
    export_fig(filename , '-eps', '-transparent', '-r900',  '-RGB');
    filename = fullfile(filepathBW, strcat(name, '.eps'));
    export_fig(filename , '-eps', '-transparent', '-r900',  '-gray');

   
end


end


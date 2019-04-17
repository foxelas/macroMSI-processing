currentLog = '2019-04-15 03_32.log';
filename = fullfile('macroMSI-processing', 'logs', '2019-04-15 03_32.log');
res = delimread(filename, ',', 'raw');
res = res.raw;

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
[~, id1] = sort(acc, 'descend');
class_s = cellfun(@(w,x,y) strrep(strjoin({w,x,y}, '|'), '_', ' ' ),...
    classifier(id1), featureSet(id1), inputSet(id1), 'uni', 0);

performance = [ auc(id1); acc(id1)];
plotClassificationPerformance('performanceComparison', 5,[], 'Performance of Top Classifiers', 'LineNames', classifier(acceptedId),'Performance', ...
    performance, 'PlotName', fullfile(savedir, 'Top_classifiers'));

performance = [auc(acceptedId); acc(acceptedId)];
plotClassificationPerformance('classificationPerformance', 6, 'LineNames', classifier(acceptedId),'Performance', ...
    performance, 'PlotName', fullfile(savedir, 'Classifiers'));
plotClassificationPerformance('classificationPerformance', 7, 'LineNames', classifier(acceptedId),'Performance', ...
    featureSet(acceptedId), 'PlotName', fullfile(savedir, 'FeatureSet'));
plotClassificationPerformance('classificationPerformance', 8,'LineNames', classifier(acceptedId),'Performance', ...
    dimred(acceptedId), 'PlotName', fullfile(savedir, 'DimRed'));
plotClassificationPerformance('classificationPerformance', 9, 'LineNames', classifier(acceptedId),'Performance', ...
    inputSet(acceptedId), 'PlotName', fullfile(savedir, 'InputSet'));
plotClassificationPerformance('classificationPerformance', 10,'LineNames', classifier(acceptedId),'Performance', ...
    {classifier(acceptedId), featureSet(acceptedId) }, 'PlotName', fullfile(savedir, 'Classifier+Features'));





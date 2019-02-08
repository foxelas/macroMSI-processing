%% 
classInput = 'unique';

classIdx = contains({classificationError.Input}, classInput);
selectedClassifiers = classificationError(classIdx);

[maxAccur, maxAccurIdx] = max([selectedClassifiers.Accuracy]);
maxAccurClassifier = selectedClassifiers(maxAccurIdx)
[maxAccurMisclassified, maxAccurFalsePositives, maxAccurFalseNegatives] = GetSelectedClassifierInfo(ID, classInput, maxAccurClassifier);

[maxAUC, maxAUCIdx] = max([selectedClassifiers.AUC]);
maxAUCClassifier = selectedClassifiers(maxAUCIdx)
[maxAUCMisclassified, maxAUCFalsePositives, maxAUCFalseNegatives] = GetSelectedClassifierInfo(ID, classInput, maxAUCClassifier);




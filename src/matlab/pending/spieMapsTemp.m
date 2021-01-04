clear all ;
load('F:\temp\mspi\matfiles\macroMSI\24-Nov-2020_Consistent_with_Veruca_Vulgaris.mat')
k = 0; 
ttestSpieMaps
load('F:\temp\mspi\matfiles\macroMSI\23-Nov-2020_Spitz_Nevus.mat')
k = 3;
ttestSpieMaps  
load('F:\temp\mspi\matfiles\macroMSI\23-Nov-2020_Melanocytic_Nevus.mat')
k = 6;
ttestSpieMaps      
load('F:\temp\mspi\matfiles\macroMSI\24-Nov-2020_Basal_Cell_Carcinoma.mat')
k = 9;
ttestSpieMaps
load('F:\temp\mspi\matfiles\macroMSI\24-Nov-2020_Pigmented_Seborrheic_Keratosis.mat')
k = 12;
ttestSpieMaps


Sample = Sample';
Method = Method';

HHbAvgUnfixedHb = HHbAvgUnfixedHb';         
HHbAvgFixedHb = HHbAvgFixedHb';
HMAvgUnfixedHb = HMAvgUnfixedHb';        
HMAvgFixedHb = HMAvgFixedHb';        
NormAvgUnfixedHb = NormAvgUnfixedHb';  
NormAvgFixedHb  =    NormAvgFixedHb';

HHbAvgUnfixedMel = HHbAvgUnfixedMel';         
HHbAvgFixedMel = HHbAvgFixedMel';
HMAvgUnfixedMel = HMAvgUnfixedMel';        
HMAvgFixedMel = HMAvgFixedMel';        
NormAvgUnfixedMel = NormAvgUnfixedMel';  
NormAvgFixedMel  =    NormAvgFixedMel';

        
% HHbStdUnfixedHb = HHbStdUnfixedHb';         
% HHbStdFixedHb = HHbStdFixedHb';
% HMStdUnfixedHb = HMStdUnfixedHb';        
% HMStdFixedHb = HMStdFixedHb';        
% NormStdUnfixedHb = NormStdUnfixedHb';  
% NormStdFixedHb  =    NormStdFixedHb';
% 
% HHbStdUnfixedMel = HHbStdUnfixedMel';         
% HHbStdFixedMel = HHbStdFixedMel';
% HMStdUnfixedMel = HMStdUnfixedMel';        
% HMStdFixedMel = HMStdFixedMel';        
% NormStdUnfixedMel = NormStdUnfixedMel';  
% NormStdFixedMel  =    NormStdFixedMel';

HMRoiMetricsHbNCC = HMRoiMetricsHbNCC';
HMRoiMetricsMelNCC = HMRoiMetricsMelNCC';
HHbRoiMetricsHbNCC = HHbRoiMetricsHbNCC'; 
HHbRoiMetricsMelNCC = HHbRoiMetricsMelNCC'; 
NormRoiMetricsHbNCC = NormRoiMetricsHbNCC'; 
NormRoiMetricsMelNCC = NormRoiMetricsMelNCC'; 
HMRoiMetricsHbHI = HMRoiMetricsHbHI';
HMRoiMetricsMelHI = HMRoiMetricsMelHI';
HHbRoiMetricsHbHI = HHbRoiMetricsHbHI'; 
HHbRoiMetricsMelHI = HHbRoiMetricsMelHI'; 
NormRoiMetricsHbHI = NormRoiMetricsHbHI'; 
NormRoiMetricsMelHI = NormRoiMetricsMelHI'; 
HMRoiMetricsHbKLD = HMRoiMetricsHbKLD';
HMRoiMetricsMelKLD = HMRoiMetricsMelKLD';
HHbRoiMetricsHbKLD = HHbRoiMetricsHbKLD'; 
HHbRoiMetricsMelKLD = HHbRoiMetricsMelKLD'; 
NormRoiMetricsHbKLD = NormRoiMetricsHbKLD'; 
NormRoiMetricsMelKLD = NormRoiMetricsMelKLD'; 
MetricsHbNCC = MetricsHbNCC';
MetricsMelNCC = MetricsMelNCC';
MetricsHbHI = MetricsHbHI';
MetricsMelHI = MetricsMelHI';
MetricsHbKLD = MetricsHbKLD';
MetricsMelKLD = MetricsMelKLD';


fullRes = table(Sample, Method, ...
HHbAvgUnfixedHb,  HHbAvgFixedHb, HMAvgUnfixedHb, HMAvgFixedHb, NormAvgUnfixedHb, NormAvgFixedHb, ...
HHbAvgUnfixedMel,  HHbAvgFixedMel, HMAvgUnfixedMel, HMAvgFixedMel, NormAvgUnfixedMel, NormAvgFixedMel, ... %HHbStdUnfixedHb,  HHbStdFixedHb, HMStdUnfixedHb, HMStdFixedHb, NormStdUnfixedHb, NormStdFixedHb, ... %HHbStdUnfixedMel,  HHbStdFixedMel, HMStdUnfixedMel, HMStdFixedMel, NormStdUnfixedMel, NormStdFixedMel, ...
     HHbRoiMetricsHbNCC, HHbRoiMetricsMelNCC, HMRoiMetricsHbNCC, HMRoiMetricsMelNCC, NormRoiMetricsHbNCC, NormRoiMetricsMelNCC, ...
        HHbRoiMetricsHbHI, HHbRoiMetricsMelHI, HMRoiMetricsHbHI, HMRoiMetricsMelHI, NormRoiMetricsHbHI, NormRoiMetricsMelHI, ...
             HHbRoiMetricsHbKLD, HHbRoiMetricsMelKLD, HMRoiMetricsHbKLD, HMRoiMetricsMelKLD, NormRoiMetricsHbKLD, NormRoiMetricsMelKLD, ...
             MetricsHbNCC, MetricsMelNCC, MetricsHbHI, MetricsMelHI, MetricsHbKLD, MetricsMelKLD);
        save('Roiresults.mat', 'fullRes');
writetable(fullRes,'Roiresults.xlsx','Sheet',1,'Range','A1')

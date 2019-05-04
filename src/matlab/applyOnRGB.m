%%Application on RGB

load('D:\temp\Google Drive\titech\research\input\saitama_v7_min_region_e\data.mat')
load('D:\temp\Google Drive\titech\research\input\saitama_v7_min_region_e\ID.mat')

dataset = 'saitama_v7_min_region_e';
skipLoading = false;
showImages = true;
saveImages = true;
tryReadData = false;

options = setOpt('ReflectanceEstimationPreset', dataset, action, showImages, saveImages, tryReadData);
options.action = 'ReflectanceEstimationPreset_rgb';
actionReflectanceEstimationComparison;
options.action = 'lbp_rgb';
actionLBP;

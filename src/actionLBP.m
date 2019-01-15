% options.systemdir = '..\..\input\saitama_v2_min_region\';
% load(fullfile(options.systemdir, 'in.mat'));
% load(fullfile(options.systemdir, 'ID.mat'));
% load(fullfile(options.systemdir, 'data.mat'));
% load(fullfile(options.systemdir, 'system.mat'));
msiN = 558; 

radius = 1;
neighbors = 8;
mapping=getmapping(neighbors,'riu2');

msibands = 7;
lbpFeatures = zeros(msiN, msibands * 10);

for k=1:msiN
	
	g = MSIStruct(k).MSI;
	mask = MSIStruct(k).Mask;
	
	gg = valueSelect(g, 'adjusted');
	
	%MATLAB built-in implementation 
	%'radius' represents scale, 'upright' represents rotation invariance
	%lbp(i) = extractLBPFeatures(gg, 'NumNeighbors', 8, 'Radius', 1, 'Upright', false)
	
    Hk = [];
    for i = 1:7
        Hk = [ Hk, lbp(squeeze(gg(i,:,:)),radius,neighbors,mapping)];
    end
    lbpFeatures(k, :) = Hk;
    
end

lbpLabels = ~[ID.IsNormal];
X = {'Benign', 'Malignant'};
lbpLabels = X(1+lbpLabels);
save('lbp.mat', 'lbpFeatures', 'lbpLabels');
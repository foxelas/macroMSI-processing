maxScale = 3;
RgbLbpFeatures = cell(maxScale, 1);
ConcatLbpFeatures = cell(maxScale, 1);
MMLbpFeatures = cell(maxScale, 1);
SumLbpFeatures = cell(maxScale, 1);

groups = findgroups([ID.MsiID]);
for g=1:max(groups)
    k = find(groups == g,1);
    files = {data([data.MsiID] == ID(k).MsiID).File}; %{data(ID(k).Data).File};   
    coordinates = [ID(groups == g).Originx; ID(groups == g).Originy];
    for lbp_operator = {'MMLBP', 'SumLBP', 'CatLBP', 'LBP'}
        lbpFeats = getLBPFeatures(lbp_operator, files, coordinates, maxScale);
    end
    if strcmp(lbp_operator, 'LBP')
        RgbLbpFeatures = addToLBPStruct(RgbLbpFeatures, lbpFeats, find(groups == g), maxScale);
    elseif strcmp(lbp_operator, 'MMLBP')
        MMLbpFeatures = addToLBPStruct(MMLbpFeatures, lbpFeats, find(groups == g), maxScale);
    elseif strcmp(lbp_operator, 'SumLBP')
        SumLbpFeatures = addToLBPStruct(SumLbpFeatures, lbpFeats, find(groups == g), maxScale);
    elseif strcmp(lbp_operator, 'CatLBP')
        ConcatLbpFeatures = addToLBPStruct(ConcatLbpFeatures, lbpFeats, find(groups == g), maxScale);
    else 
        disp('Unrecognized LBP method. Aborting...')
        return      
    end

end

save( fullfile(options.saveOptions.savedir, 'ReflectanceEstimationPreset', 'out.mat'),'ConcatLbpFeatures', 'SumLbpFeatures', 'MMLbpFeatures','RgbLbpFeatures', '-append');

    
function [current] = addToLBPStruct(base, lbps, ks, scales)
    for sc=1:scales
        for k=1:length(ks)
            base{sc,1}(ks(k), :) = lbps{sc}(k,:); 
        end
    end
    current = base;
end


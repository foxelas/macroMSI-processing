maxScale = 3;
RgbLbpFeatures = cell(maxScale, 1);
ConcatLbpFeatures = cell(maxScale, 1);
MMLbpFeatures = cell(maxScale, 1);
SumLbpFeatures = cell(maxScale, 1);

for k = 1:msiN
    infile = fullfile(options.systemdir, 'infiles', strcat('poi_', num2str(k), '.mat'));
    load(infile, 'poiName', 'poiRAW', 'poiSegmentMask', ...
        'roiSeeds', 'measuredSpectrum', 'poiWhite');
    for lbp_operator = {'MMLBP', 'SumLBP', 'CatLBP', 'LBP'}
        lbpFeats = getLBPFeatures(lbp_operator, options, k, maxScale);

        if strcmp(lbp_operator, 'LBP')
            RgbLbpFeatures = addToLBPStruct(RgbLbpFeatures, lbpFeats, k, maxScale);
        elseif strcmp(lbp_operator, 'MMLBP')
            MMLbpFeatures = addToLBPStruct(MMLbpFeatures, lbpFeats, k, maxScale);
        elseif strcmp(lbp_operator, 'SumLBP')
            SumLbpFeatures = addToLBPStruct(SumLbpFeatures, lbpFeats, k, maxScale);
        elseif strcmp(lbp_operator, 'CatLBP')
            ConcatLbpFeatures = addToLBPStruct(ConcatLbpFeatures, lbpFeats, k, maxScale);
        else
            disp('Unrecognized LBP method. Aborting...')
            return
        end
    end

end

filename = mkNewDir(fullfile(options.saveOptions.savedir, getOutputDirectoryMap('features'), 'out.mat'));
if exist(filename, 'file')
    save(filename, 'ConcatLbpFeatures', 'SumLbpFeatures', 'MMLbpFeatures', 'RgbLbpFeatures', '-append');
else
    save(filename, 'ConcatLbpFeatures', 'SumLbpFeatures', 'MMLbpFeatures', 'RgbLbpFeatures');
end

if (options.showImages)
    for g = 1:max([ID.Group])
        infile = fullfile(options.systemdir, 'infiles', strcat('group_', num2str(g), '.mat'));
        load(infile, 'poiName', 'raw', 'whiteReference', 'specimenMask');
        visualizeLBP(raw, whiteReference, specimenMask, g, options.saveOptions);
    end
end

function [current] = addToLBPStruct(base, lbps, ks, scales)
for sc = 1:scales
    for k = 1:length(ks)
        base{sc, 1}(ks(k), :) = lbps{sc}(k, :);
    end
end
current = base;
end

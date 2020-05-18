function [lbpFeats] = getLBPFeatures(type, options, k, maxScale, neighbors, mapping)
%main code provided by http://www.cse.oulu.fi/CMV/Downloads/LBPMatlab
    if (nargin < 4)
        maxScale = 3;
    end
    if (nargin < 5)
        neighbors = 8;
    end
    if (nargin < 6)
        mapping=getmapping(neighbors,'riu2');
    end

    msibands = 9;
    %spectralNeighbors = 2;
    ConcatLbpFeatures = cell(maxScale, 1);
    SumLbpFeatures = cell(maxScale, 1);
    MMLbpFeatures = cell(maxScale, 1);
    RgbLbpFeatures = cell(maxScale, 1);
    
    for scale = 1:maxScale
        riubins = (neighbors + 2);
        concatlbpFeatures = zeros(1, msibands * riubins);
        sumlbpFeatures = zeros(1, riubins);
        mmlbpFeatures = zeros(1, 4 * riubins);
        rgblbpFeatures = zeros(1, riubins);

        infile = fullfile(options.systemdir, 'infiles', strcat('poi_', num2str(k), '.mat'));
        load(infile,  'poiRAW','poiWhite');
        msi = raw2msi(poiRAW, 'extended');
        rgb = poiWhite;

        if strcmp(type, 'CatLBP') || strcmp(type, 'SumLBP')
            for i = 1:msibands
                lbps = lbp(squeeze(msi(i,:,:)),scale,neighbors,mapping);
                if strcmp(type, 'CatLBP'); concatlbpFeatures((i-1)*10 + (1:10)) = lbps; end
                if strcmp(type, 'SumLBP'); sumlbpFeatures = sumlbpFeatures + lbps; end
            end
        elseif strcmp(type, 'MMLBP')
            mmlbpFeatures = lbp(msi, scale, neighbors, mapping);
        elseif strcmp(type, 'LBP')
            gr = rgb2gray(rgb);
            rgblbpFeatures = lbp(gr,scale,neighbors,mapping);
        end
        
        ConcatLbpFeatures{scale} = concatlbpFeatures;
        SumLbpFeatures{scale} = sumlbpFeatures;
        MMLbpFeatures{scale} = mmlbpFeatures;
        RgbLbpFeatures{scale} = rgblbpFeatures ;        
    end
    
    if strcmp(type, 'CatLBP')
        lbpFeats = ConcatLbpFeatures;
    elseif strcmp(type, 'SumLBP')
        lbpFeats = SumLbpFeatures;
    elseif strcmp(type, 'MMLBP')
        lbpFeats = MMLbpFeatures;
    else 
        lbpFeats = RgbLbpFeatures;
    end
end
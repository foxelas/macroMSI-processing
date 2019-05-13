function [lbpFeats] = getLBPFeatures(type, files, coordinates, maxScale, neighbors, mapping)
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
    spectralNeighbors = 2;
    ConcatLbpFeatures = cell(maxScale, 1);
    SumLbpFeatures = cell(maxScale, 1);
    MMLbpFeatures = cell(maxScale, 1);
    RgbLbpFeatures = cell(maxScale, 1);

    if ~(size(coordinates,2) == 2);  coordinates = coordinates'; end
    rois = size(coordinates, 1);
    for scale = 1:maxScale
        riubins = (neighbors + 2);
        concatlbpFeatures = zeros(rois, msibands * riubins);
        sumlbpFeatures = zeros(rois, riubins);
        mmlbpFeatures = zeros(rois, 4 * riubins);
        rgblbpFeatures = zeros(rois, riubins);

        for k = 1:rois
            coor = coordinates(k,:);

            [msi, rgb] = readMSI(files, coor, 5 + scale, 5 + scale, []); 
            msi = raw2msi(msi, 'extended');

            if strcmp(type, 'CatLBP') || strcmp(type, 'SumLBP')
                for i = 1:msibands
                    lbps = lbp(squeeze(msi(i,:,:)),scale,neighbors,mapping);
                    if strcmp(type, 'CatLBP'); concatlbpFeatures(k, (i-1)*10 + (1:10)) = lbps; end
                    if strcmp(type, 'SumLBP'); sumlbpFeatures(k,:) = sumlbpFeatures(1,:) + lbps; end
                end
            elseif strcmp(type, 'MMLBP')
                mmlbpFeatures(k,:) = lbp(msi, scale, neighbors, mapping);
            elseif strcmp(type, 'LBP')
                gr = rgb2gray(rgb);
                rgblbpFeatures(k,:) = lbp(gr,scale,neighbors,mapping);
            end
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
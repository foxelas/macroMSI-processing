%code provided by http://www.cse.oulu.fi/CMV/Downloads/LBPMatlab
    
maxScale = 3;
neighbors = 8;
spectralNeighbors = 2;
mapping=getmapping(neighbors,'riu2');

if ~contains(options.action, 'rgb') 
    msibands = 9;
    ConcatLbpFeatures = cell(maxScale, 1);
    SumLbpFeatures = cell(maxScale, 1);
    MMLbpFeatures = cell(maxScale, 1);

    for scale = 1:maxScale
        riubins = (neighbors + 2);
        concatlbpFeatures = zeros(msiN, msibands * riubins);
        sumlbpFeatures = zeros(msiN, riubins);
        mmlbpFeatures = zeros(msiN, 4 * riubins);

        for k=1:msiN

            files = {data([data.MsiID] == ID(k).MsiID).File}; %{data(ID(k).Data).File};   
            coordinates = [ID(k).Originx; ID(k).Originy];

            msi = readMSI(files, coordinates, 5 + scale, 5 + scale, []); 
            gg = raw2msi(msi, 'extended');

            for i = 1:msibands
                lbps = lbp(squeeze(gg(i,:,:)),scale,neighbors,mapping);
                concatlbpFeatures(k, (i-1)*10 + (1:10)) = lbps;
                sumlbpFeatures(k,:) = sumlbpFeatures(k,:) + lbps;
            end

            mmlbpFeatures(k,:) = lbp(gg, scale, neighbors, mapping);
        end

        ConcatLbpFeatures{scale} = concatlbpFeatures / max(concatlbpFeatures(:));
        SumLbpFeatures{scale} = sumlbpFeatures / max(sumlbpFeatures(:));
        MMLbpFeatures{scale} = mmlbpFeatures / max(mmlbpFeatures(:));
    end
    save( fullfile(options.saveOptions.savedir, 'ReflectanceEstimationPreset', 'out.mat'),'ConcatLbpFeatures', 'SumLbpFeatures', 'MMLbpFeatures', '-append');
else
     RgbLbpFeatures = cell(maxScale, 1);
    for scale = 1:maxScale
        riubins = (neighbors + 2);
        rgblbpFeatures = zeros(msiN, riubins);
        for k=1:msiN

            files = {data([data.MsiID] == ID(k).MsiID).File}; %{data(ID(k).Data).File};   
            coordinates = [ID(k).Originx; ID(k).Originy];

            [~, whiteI] = readMSI(files, coordinates, 5 + scale, 5 + scale, []); 
            whiteI = rgb2gray(whiteI);
            rgblbpFeatures(k,:) = lbp(whiteI,scale,neighbors,mapping);
        end

        RgbLbpFeatures{scale} = rgblbpFeatures / max(rgblbpFeatures(:));
    end
    save( fullfile(options.saveOptions.savedir, 'ReflectanceEstimationPreset', 'out.mat'),'RgbLbpFeatures', '-append');
end

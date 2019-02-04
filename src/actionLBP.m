%code provided by http://www.cse.oulu.fi/CMV/Downloads/LBPMatlab
out = matfile( fullfile(options.saveOptions.savedir, 'ReflectanceEstimationPreset', 'out.mat'), 'Writable', true);

maxScale = 3;
neighbors = 8;
mapping=getmapping(neighbors,'riu2');

msibands = 7;
multiScaleLbpFeatures = cell(maxScale, 1);
for scale = 1:maxScale
    
    lbpFeatures = zeros(msiN, msibands * 10);
    for k=1:msiN
        
        files = {data(ID(k).Data).File};   
        coordinates = [ID(k).Originx; ID(k).Originy];
        % coordinates = out.newCoordinates(k,:);
        segment = readMSI(files, coordinates, 5 + scale, 5 + scale, []); 

        g = segment.MSI;
        gg = raw2msi(g, 'adjusted');

        for i = 1:7
            lbpFeatures(k, (i-1)*10 + (1:10)) = lbp(squeeze(gg(i,:,:)),scale,neighbors,mapping);
        end
    end
    
    multiScaleLbpFeatures{scale} = lbpFeatures / max(lbpFeatures(:));

end
out.multiScaleLbpFeatures = multiScaleLbpFeatures;

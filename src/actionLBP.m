%code provided by http://www.cse.oulu.fi/CMV/Downloads/LBPMatlab

maxScale = 2;
neighbors = 8;
mapping=getmapping(neighbors,'riu2');

msibands = 7;
MultiScaleLbpFeatures = cell(maxScale, 1);
for scale = 1:maxScale
    
    lbpFeatures = zeros(msiN, msibands * 10);
    for k=1:msiN
        
        files = {data(ID(k).Data).File};   
        coordinates = [ID(k).Originx; ID(k).Originy];
        % coordinates = out.newCoordinates(k,:);
        msi = readMSI(files, coordinates, 5 + scale, 5 + scale, []); 
        gg = raw2msi(msi, 'adjusted');

        for i = 1:7
            lbpFeatures(k, (i-1)*10 + (1:10)) = lbp(squeeze(gg(i,:,:)),scale,neighbors,mapping);
        end
    end
    
    MultiScaleLbpFeatures{scale} = lbpFeatures / max(lbpFeatures(:));

end
save( fullfile(options.saveOptions.savedir, 'ReflectanceEstimationPreset', 'out.mat'),'MultiScaleLbpFeatures', '-append');
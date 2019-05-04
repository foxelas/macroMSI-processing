        %% Reconstruct sRGB from MSI
        
        for k = [81]
            [~, sampleName] = generateName(options, 'plot', ID(k));
            [g , whiteImg ]= readMSI({data([data.MsiID] == ID(k).MsiID).File});
            figure(5); imshow(whiteImg);
            sRGB = createSRGB(g, 'original', sampleName, ID(k), options);
        end
        
        % end of create srgb
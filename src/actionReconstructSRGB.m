        %% CREATE SRGB
        
        for k = [250, 270]
            [~, sampleName] = generateName(options, 'plot+save', data(ID(k).Representative), ID(k));
            g = readMSI({data(ID(k).Data).File}, [], [], [], [], [0, 450, 465, 505, 525, 575, 605, 630], true);
            sRGB = createSRGB(g, 'original', sampleName, ID(k), options);
        end
        
        % end of create srgb
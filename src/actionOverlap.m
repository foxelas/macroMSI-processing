        %% Overlap of cancer / fixed
        
        if containts(options.action, 'Measured')
            [~, xx] = unique(strcat({ID.Csvid}, {ID.T}), 'last');
            id = ID(xx);
            curves = zeros(length(380:780), length(id));
        else
            e = matfile(fullfile(options.saveOptions.saveDir, 'ReflectanceEstimationSimple', 'out.mat'));
            id = ID;
            curves = e.EstimatedSpectrumStruct;
        end
        
        %cancer
        lineNamesCancer = cell(length(id), 1);
        lineNamesFixed = cell(length(id), 1);
        lineNamesSample = cell(length(id), 1);
        
        for k = 1:length(id)
            curves(:, k) = measuredSpectrumStruct(k).Spectrum;
            if id(k).IsNormal
                lineNamesCancer{k} = 'normal';
            else
                lineNamesCancer{k} = 'cancerous';
            end
            if id(k).IsFixed
                lineNamesFixed{k} = 'fixed';
            else
                lineNamesFixed{k} = 'unfixed';
            end
            if id(k).IsNormal
                lineNamesSample{k} = [id(k).Sample, '_normal'];
            else
                lineNamesSample{k} = [id(k).Sample, '_cancerous'];
            end
        end
        options.saveOptions.plotName = generateName(options, [options.action, 'Cancer']);
        plots('overlapSpectrum', 2, 380:780, curves, 'Cancer', 'wavelength', wavelength, 'lineNames', lineNamesCancer, 'saveOptions', options.saveOptions);
        options.saveOptions.plotName = generateName(options, [options.action, 'Fixed']);
        plots('overlapSpectrum', 3, 380:780, curves, 'Fixed', 'wavelength', wavelength, 'lineNames', lineNamesFixed, 'saveOptions', options.saveOptions);
        options.saveOptions.plotName = generateName(options, [options.action, 'Sample']);
        plots('overlapSpectrum', 1, 380:780, curves, 'Sample', 'wavelength', wavelength, 'lineNames', lineNamesSample, 'saveOptions', options.saveOptions);
        
        % end overlap of cancer / fixed
width = 3;
height = 3;
ID2 = load('..\MATLAB\Data\saitamav2\ID_by_hand.mat');

showRef = false;
showImage = false;

for k = [66, 165, 175, 236, 297]
    if ~(ID(k).IsFixed) && contains(data(ID(k).Representative).Type, 'unfixed')
        
        %% Choose search step
        info = imfinfo(data(ID(k).Representative).File);
        imheight = info.Height;
        imwidth = info.Width;
        searchboxR = min([ceil(max(imwidth, imheight)/70), 20]);
        fprintf('Now optimization at sample ID(%d), searching inside square with edge %d.\n', k, 2*searchboxR);
        xx = ID(k).Originx;
        yy = ID(k).Originy;
        g = readMSI({data(ID(k).Data).File}, xx, yy, width, height);
        
        % Measured spectrum
        measured = interp1(380:780, eval(strcat('in.s_', generateName([], 'csv', [], ID(k)))), wavelength, 'nearest');
        
        %% Error of the hand selected (old) point
        % Retrieve correction coefficients
        %          options.coeff = squeeze(coeff(k, min(find(strcmp(pixelValueSelectionMethods , options.pixelValueSelectionMethod)),3), :))';
        
        [~, rold, rmseold, ~] = reflectanceEstimation(g, measured, ID(k), options);
        rmsemin = rmseold;
        
        if (showImage)
            showCroppedSection([], data(ID(k).Representative).File, xx, yy, [], [], true, 'y');
        end
        
        measuredBatch = measured;
        
        %% Create square searchbox
        for x = (xx - searchboxR):(xx + searchboxR)
            for y = (yy - searchboxR):(yy + searchboxR)
                
                % read new image and do ref est
                g = readMSI({data(ID(k).Data).File}, x, y, width, height);
                [~, rnew, rmsenew, ~] = reflectanceEstimation(g, measured, ID(k), options);
                
                % compare with previous error
                if (rmsenew < rmsemin)
                    xxmin = x;
                    yymin = y;
                    rmsemin = rmsenew;
                    rmin = rnew;
                end
                
                measuredBatch = [measuredBatch, rnew];
                
                %visual data
                if (showRef)
                    plots('estimationComparison', 2, [measured, rnew, rold], '', 'wavelength', wavelength, 'lineNames', {'MS center \lambda', 'Measured', 'estimated(current error)', 'estimated(original)'}, 'saveOptions', options.saveOptions);
                    plots('singlemeasurement', 3, rnew, '', 'wavelength', wavelength, 'lineNames', 'Candidate sample', 'saveOptions', options.saveOptions);
                end
                if (showImage)
                    showCroppedSection([], data(ID(k).Representative).File, x, y, [], [], false, 'm');
                end
                
            end
        end
        
        %% Check if accept the new min or not
        if (abs(ID(k).Originx-xxmin) < 30 && abs(ID(k).Originy-yymin) < 30)
            ID(k).Originx = xxmin;
            ID(k).Originy = yymin;
            ID(k).rmse = rmsemin;
        end
        
        ID(k).IsFixed = true;
        idmin(k, 1:2) = [xxmin, yymin];
        
        %% To check before and after estimation
        hold on
        plots('allEstimations', 2, [measuredBatch, rmin]', '', 'wavelength', wavelength, 'saveOptions', options.saveOptions);
        hold off
        title(sprintf('Original error %.3f, Minimum error %.3f, difference %.3f . Discard selection by-hand: %d', rmseold, rmsemin, abs(rmsemin-rmseold), rmsemin ~= rmseold))
        
        
        showCroppedSection([], data(ID(k).Representative).File, xx, yy, [], [], true, 'y');
        showCroppedSection([], data(ID(k).Representative).File, xxmin, yymin, [], [], false, 'm');
        
        
        %         hold on
        %         gg = valueSelect(gNew, 'rms');
        %         ref = mean(mean(gg,3),2);
        %         plot(bands(2:end), ref', 'm*-'); %magenda new
        %         gg = valueSelect(gOld, 'rms');
        %         ref = mean(mean(gg,3),2);
        %         plot(bands(2:end), ref', 'c*-'); %cyan old
        %         hold off
        %         pause(0.05)
    end
    
    save('IDmin.mat', 'ID');
    
end

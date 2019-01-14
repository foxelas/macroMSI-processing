systemdir = {'D:\temp\Google Drive\titech\research\input\saitama_v2_min_square\', 'D:\temp\Google Drive\titech\research\input\saitama_v2_min_region\'};
for i = 2:2
    
    options.systemdir = systemdir{i}
    load(fullfile(options.systemdir, 'in.mat'));
    load(fullfile(options.systemdir, 'ID.mat'));
    load(fullfile(options.systemdir, 'data.mat'));

    pixelValueSelectionMethods = {'green', 'rms', 'adjusted'};

    % types = unique({ID.Type});
    % samples = unique({ID.Sample});
    % for type = types
    %     for sample = samples
    %         points = unique( {ID(find(strcmp({ID.Sample}, sample) & strcmp({ID.Type}, type) )).IMG} );

    %c = matfile([options.systemdir, 'coeff.mat'], 'Writable', true);
    coeff = ones(length(ID), 3, 7);
    src = '../../input/others/coeff.xlsx';
    for k = 202:length(ID)
        sampleName = generateName([], 'image', data(ID(k).Representative), ID(k));

        % Retrieve MSI data
        g = MSIStruct(k).MSI;

        % Retrieve spectrum data
        refmeasured = MeasuredSpectrumStruct(k).Spectrum;

        for m = 1:length(pixelValueSelectionMethods)
            gg = valueSelect(g, pixelValueSelectionMethods{m});
            if contains(options.systemdir, 'region')
                mask = reshape( MSIStruct(k).Mask, 1, size(gg,2)*size(gg,3));
                gg = reshape(gg, size(gg,1), size(gg,2)*size(gg,3));
                gg = gg(:,mask);
            else
                gg = gg(:,2:4,2:4); %optional step
                gg = reshape(gg, size(gg,1), 3*3);
            end
            refmsi = im2uint16(mean(gg, 2))';
            xlswrite(src, refmeasured, 4, 'B3');
            xlswrite(src, refmsi, 4, 'D407:J407');
            ctemp = xlsread(src, 4, 'D414:J414');
            coeff(k, m, :) = ctemp;
        end
        k
    end

    save([options.systemdir, 'coeff.mat'], 'coeff', '-v7.3');
end


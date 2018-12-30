pixelValueSelectionMethods = {'green', 'rms', 'adjusted'};

% types = unique({ID.Type});
% samples = unique({ID.Sample});
% for type = types
%     for sample = samples
%         points = unique( {ID(find(strcmp({ID.Sample}, sample) & strcmp({ID.Type}, type) )).IMG} );

%c = matfile([options.systemdir, 'coeff.mat'], 'Writable', true);
coeff = ones(length(ID), 3, 7);

for k = 1:length(ID)
    sampleName = generateName([], 'image', data(ID(k).Representative), ID(k));
    
    % Retrieve MSI data
    g = MSIStruct(k).MSI;
    
    % Retrieve spectrum data
    refmeasured = measuredSpectrumStruct(k).Spectrum;
    
    for m = 1:length(pixelValueSelectionMethods)
        gg = valueSelect(g, pixelValueSelectionMethods{m});
        refmsi = im2uint16(mean(mean(gg, 3), 2))';
        xlswrite('coeff.xlsx', refmeasured, 4, 'B3');
        xlswrite('coeff.xlsx', refmsi, 4, 'D407:J407');
        ctemp = xlsread('coeff.xlsx', 4, 'D414:J414');
        coeff(k, m, :) = ctemp;
    end
    k
end

save([options.systemdir, 'coeff.mat'], 'coeff');
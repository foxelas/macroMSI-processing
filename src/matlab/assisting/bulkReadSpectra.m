function [Spectra,CompleteSpectra,SpectraNames] = bulkReadSpectra(ID, data, wavelengthN, options)
    %% Read raw spectra 
    fprintf('Reading spectral data according to ID file.\n');

    msiN = length(ID);

    wavelengthStep = ceil(length(380:780) / wavelengthN);
    wavelengthIdxs = 1:wavelengthStep:length(380:780);
    [~, uniqueSpectraIdxs, uniqueSpectraIdxsInID] = unique(strcat({ID.SpectrumFile}, {ID.T}));
    idd = ID(uniqueSpectraIdxs);
    specN = length(uniqueSpectraIdxs);
    completeUniqueSpectra = zeros(specN, 401);

    for i = 1:specN
        % read raw measured spectrum
        rawSpectrum = readSpectrum(idd(i).SpectrumFile, idd(i).T);
        % read raw white measured spectrum of the reference surface
        referenceSpectrum = readSpectrum(char(strcat(data(idd(i).RgbID).Sample, '\', 'white.csv')));

        if abs(rawSpectrum-referenceSpectrum) < 0.000001
            error('Measurement is same as white.')
        end
        completeUniqueSpectra(i,:) = rawSpectrum ./ referenceSpectrum;  
    end
    uniqueSpectra = completeUniqueSpectra(:,wavelengthIdxs);
    clear('idd');

    %% Save Spectra Struct
    Spectra = zeros(msiN, wavelengthN);
    CompleteSpectra = zeros(msiN, 401);
    SpectraNames = cell(msiN, 1);
    for i = 1:msiN               
        Spectra(i,:) = uniqueSpectra(uniqueSpectraIdxsInID(i),:);
        CompleteSpectra(i,:) = completeUniqueSpectra(uniqueSpectraIdxsInID(i),:);
        SpectraNames{i} = strcat(strrep(ID(i).SpectrumFile, '\', '_'), ', ', ID(i).T);
    end 
    save(generateName('matfilein', options), 'Spectra', 'CompleteSpectra', ...
        'SpectraNames', '-append');

end


function [ ] = readAndSaveSpectrum( options, ID, data )
%READANDSAVESPECTRUM   Reads and saves the normalized spectrum for all samples in ID
% savedir: the source of the data to be read
% ID,data: functions containing information about the data to me read 

fprintf('Reading spectral data according to ID file.\n');

measuredSpectrumStruct = struct('Name', {}, 'Index', [], 'Spectrum', []);
appendingIdx = 1;

if ~isempty(strfind(ID(1).Csvid,'.csv'))
    % has csv measured info 
    for i  = 1:length(ID)
        
        % read raw measured spectrum
        measuredSpectrum = readSpectrum(ID(i).Csvid, ID(i).T);
        
        % read raw white measured spectrum of the reference surface 
        white = readSpectrum( char(strcat(data(ID(i).Representative).Sample, '\', 'white.csv') ) );
        
        if abs(measuredSpectrum  - white) < 0.000001
            i
            error('Measurement is same as white.')
        end
        
        if  ~(isempty(measuredSpectrum))    
            measuredSpectrumStruct(appendingIdx) = struct('Name', generateName([], 'csvfile', [], ID(i)), 'Index', appendingIdx, 'Spectrum', measuredSpectrum ./ white);
            appendingIdx = length(measuredSpectrumStruct) + 1;
        end

    end
    
elseif ~isempty(dir([options.systemdir, '*MeasuredReflectance.mat']))
    % has mat measured info 
    measured = dir([options.systemdir, '*MeasuredReflectance.mat']);
    load(measured.name);
    for i  = 1:length(ID)
        idx1 = ID(i).Representative;
        idx2 = find(strcmp({measuredReflectance.sampleName}, data(idx1).Sample));
        measuredSpectrumStruct(appendingIdx) = struct('Name', generateName([], 'csvfile', [], ID(i)), 'Index', appendingIdx, 'Spectrum', measuredReflectance(idx2).reflectance);
        appendingIdx = length(measuredSpectrumStruct) + 1;

    end
end

m = matfile(generateName(options, 'matfilein'),'Writable',true);
m.MeasuredSpectrumStruct = measuredSpectrumStruct;

fprintf('Finished reading all spectral data according to ID file. Saving in %s.\n', options.systemdir);
    
end


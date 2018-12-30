function [] = readAndSaveMSI(options, ID, data, width, height, bands)

%% READANDSAVEMSI reads all files in ID and saves them in .mat file
% datadir: directory to fetch data to be read
% savedir: directory to save data after read

%    % fieldnames(ID);
%     str = sprintf('%f,', ID.originx);
%     originx = textscan(str, '%f', 'delimiter', ',', 'EmptyValue', NaN);
%     ID = ID(~isnan(originx{:})); %All MSI samples have 1) spectrogram point position 2) measured spectrum
%     save(strcat(savedir, '\ID.mat'), 'ID');
fprintf('Reading MSI data according to ID file.\n');

if (nargin < 4)
    width = 5;
end

if (nargin < 5)
    height = 5;
end

if (nargin < 6)
    bands = [0, 450, 465, 505, 525, 575, 605, 630];
end

whiteMSIStruct = struct('Name', {}, 'Index', [], 'MSI', []);
darkMSIStruct = struct('Name', {}, 'Index', [], 'MSI', []);
MSIStruct = struct('Name', {}, 'Index', [], 'MSI', []);
appendingIdx = 1;

for i = 1:length(ID)
    
    [MSI, whiteReference, darkReference] = readMSI({data(ID(i).Data).File}, ID(i).Originx, ID(i).Originy, width, height, bands, '..\MATLAB\Cropped\', [ID(i).Csvid, num2str(ID(i).UniqueCount)]);
    
    name = generateName([], 'image', data(ID(i).Representative), ID(i));
    %             name = strjoin( [ data(idx(1)).Camera, strrep( strrep(ID(i).csvid, '.csv', ''), '\', '_') , num2str(ID(i).UniqueCount)] ,'_');
    
    if ~(isempty(whiteReference))
        whiteMSIStruct(appendingIdx) = struct('Name', name, 'Index', appendingIdx, 'MSI', whiteReference);
    end
    
    if ~(isempty(darkReference))
        darkMSIStruct(appendingIdx) = struct('Name', name, 'Index', appendingIdx, 'MSI', darkReference);
    end
    
    if ~(isempty(MSI))
        MSIStruct(appendingIdx) = struct('Name', name, 'Index', appendingIdx, 'MSI', MSI);
        appendingIdx = length(MSIStruct) + 1;
    end
end

m = matfile(generateName(options, 'matfilein'), 'Writable', true);
m.MSIStruct = MSIStruct;
m.WhiteMSIStruct = whiteMSIStruct;
m.DarkMSIStruct = darkMSIStruct;

fprintf('Finished reading all MSI data according to ID file. Saving in %s.\n', options.systemdir);

end

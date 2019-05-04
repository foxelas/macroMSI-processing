function [dirCurrentName, currentName] = generateName(options, filename, idx)
%GENERATEOUTNAME generates the name of the file,
%where the result will be saved.
%savefilename: the saveas filename icluding directory in the form (..\..\..\savename.extension)
%sampleName: the sample name

if (nargin < 3)
    idx = [];
end
    
currentName = [];
dirCurrentName = [];

switch filename
    case 'matfilein'
        dirCurrentName = fullfile(options.systemdir, 'in.mat');
        
    case 'matfilein-v7.3'
        dirCurrentName = fullfile(options.systemdir, 'in-v73.mat');

    case 'matfileout'
        dirCurrentName = fullfile(options.saveOptions.savedir, options.action, 'out.mat');
        directory = fileparts(dirCurrentName);
        mkdir_custom(directory);
        
    case 'action detail'
        [~, currentName] = generateName(options, 'current', idx);

        extra = {options.smoothingMatrixMethod, options.noiseType};
        if strcmp(options.action, 'ReflectanceEstimation')
            if strcmp(options.smoothingMatrixMethod, 'markovian')
                extra{end+1} = num2str(options.rho);
            end
        end
        extra{end+1} = currentName;
        currentName = strjoin(extra, '_');
        dirCurrentName = fullfile(options.saveOptions.savedir, options.action, currentName);

    case 'csv'
        if isempty(idx.SpectrumFile)
            currentName = num2str(idx.MsiID);
        else
            currentName = strrep(idx.SpectrumFile, '.csv', '');
            currentName = strrep(currentName, '\', '_');
        end

    case 'sample'
        splits = strsplit(idx.SpectrumFile, '\');
        currentName = splits{1};

    case 'current'  
        [~, csv] = generateName(options, 'csv', idx);
        time = strrep(idx.T, ':', ' ');
        time = strrep(time, '.', ',');
        un = strcat('(', num2str(idx.MsiID) , ')');
        load(fullfile(options.systemdir, 'data.mat'), 'data');
        datax = data(:, idx.RgbID);
        if (~isempty(datax) && isnumeric(datax.Sample))
            datax.Sample = num2str(datax.Sample);
        end
        currentName = strjoin({csv, time, char(datax.Camera), un }, '_');

    case 'read'    
        [~, currentName] = generateName(options, 'current', idx);
        dirCurrentName = fullfile(options.saveOptions.savedir, 'Cropped', currentName);
        directory = fileparts(dirCurrentName);
        mkdir_custome(directory);

    case 'plot'
        [~, currentName] = generateName(options, 'current', idx);
        dirCurrentName = fullfile(options.saveOptions.savedir, options.action, currentName);
        directory = fileparts(dirCurrentName);
        mkdir_custom(directory);

    case 'input'
        dirCurrentName = fullfile('..', '..', '..', 'input');
        
    case 'output'
        dirCurrentName = fullfile('..', '..', '..', 'output');
        
    otherwise
        dirCurrentName = fullfile(options.saveOptions.savedir, options.action, filename);
        directory = fileparts(dirCurrentName);
        mkdir_custome(directory);

end

end

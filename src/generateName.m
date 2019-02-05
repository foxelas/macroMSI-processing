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
        if ~exist(directory, 'dir')
            mkdir(directory);
            addpath(directory);
        end
        
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
        if isempty(idx.Csvid)
            currentName = num2str(idx.UniqueCount);
        else
            currentName = strrep(idx.Csvid, '.csv', '');
            currentName = strrep(currentName, '\', '_');
        end

    case 'sample'
        splits = strsplit(idx.Csvid, '\');
        currentName = splits{1};

    case 'current'  
        [~, csv] = generateName(options, 'csv', idx);
        time = strrep(idx.T, ':', ' ');
        time = strrep(time, '.', ',');
        un = strcat('(', num2str(idx.UniqueCount) , ')');
        load(fullfile(options.systemdir, 'data.mat'), 'data');
        datax = data(:, idx.Representative);
        if (~isempty(datax) && isnumeric(datax.Sample))
            datax.Sample = num2str(datax.Sample);
        end
        currentName = strjoin({csv, time, char(datax.Camera), un }, '_');

    case 'read'    
        [~, currentName] = generateName(options, 'current', idx);
        dirCurrentName = fullfile(options.saveOptions.savedir, 'Cropped', currentName);
        directory = fileparts(dirCurrentName);
        if ~exist(directory, 'dir')
            mkdir(directory);
            addpath(directory);
        end

    case 'plot'
        [~, currentName] = generateName(options, 'current', idx);
        dirCurrentName = fullfile(options.saveOptions.savedir, options.action, currentName);
        directory = fileparts(dirCurrentName);
        if ~exist(directory, 'dir')
            mkdir(directory);
            addpath(directory);
        end

    otherwise
        dirCurrentName = fullfile(options.saveOptions.savedir, options.action, filename);
        directory = fileparts(dirCurrentName);
        if ~exist(directory, 'dir')
            mkdir(directory);
            addpath(directory);
        end
end

end

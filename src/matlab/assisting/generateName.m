function name = generateName(namingCase, options, idx)
%GENERATEOUTNAME generates the name of the file,
%where the result will be saved.
% Inputs:
% namingCase - the naming case or specific filename 
% options - running options
% idx - the information of present data sample
% 
% Outputs: 
% dirCurrentName - directory only 
% 
if (nargin < 3)
    idx = [];
end
    
name = [];

switch namingCase
    case 'matfilein'
        name = fullfile(options.systemdir, 'in.mat');
        
    case 'matfilein-v7.3'
        name = fullfile(options.systemdir, 'in-v73.mat');

    case 'matfileout'
        name = fullfile(options.saveOptions.savedir, options.action, 'out.mat');
        
%     case 'action detail'
%         name = generateName(options, 'current', idx);
% 
%         extra = {options.smoothingMatrixMethod, options.noiseType};
%         if strcmp(options.action, 'ReflectanceEstimation')
%             if strcmp(options.smoothingMatrixMethod, 'markovian')
%                 extra{end+1} = num2str(options.rho);
%             end
%         end
%         extra{end+1} = name;
%         name = strjoin(extra, '_');
%         name = fullfile(options.saveOptions.savedir, options.action, name);

    case 'csv'
        if isempty(idx.SpectrumFile)
            name = num2str(idx.MsiID);
        else
            name = strrep(idx.SpectrumFile, '.csv', '');
            name = strrep(name, '\', '_');
        end

    case 'sample'
        splits = strsplit(idx.SpectrumFile, '\');
        name = splits{1};

    case 'current'  
        csv = generateName('csv', options, idx);
        time = strrep(idx.T, ':', ' ');
        time = strrep(time, '.', ',');
        un = strcat('(', num2str(idx.MsiID) , ')');
        sample = idx.Sample;
%         load(fullfile(options.systemdir, 'data.mat'), 'data');
%         datax = data(:, idx.RgbID);
%         if (~isempty(datax) && isnumeric(datax.Sample))
%             datax.Sample = num2str(datax.Sample);
%         end
        name = strjoin({csv, time, sample, un }, '_');

    case 'read'    
        currentName = generateName('current', options, idx);
        name = fullfile(options.saveOptions.savedir, 'Cropped', currentName);

    case 'plot'
        currentName = generateName('current', options, idx);
        name = fullfile(options.saveOptions.savedir, options.action, currentName);

    case 'input'
        if isempty(options)
            name = fullfile('..', '..', '..', 'input');
        else
            name = options.systemdir;
        end
        
    case 'output'
        if isempty(options)
            name = fullfile('..', '..', '..', 'output');
        else
            name = fullfile(options.saveOptions.savedir, options.action);
        end
        
    case 'id'
        name = fullfile(options.systemdir, 'ID.mat');
        
    case 'system'
        name = fullfile(options.systemdir, 'system.mat');
        
    otherwise
        name = fullfile(options.saveOptions.savedir, options.action, namingCase);

end

directory = fileparts(name);
mkdir_custom(directory);

end

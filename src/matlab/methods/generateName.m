function name = generateName(namingCase, idx)
%   GENERATENAME generates the name of the file,
%   where the result will be saved.
%   Inputs:
%   namingCase - the naming case or specific filename
%   idx - the information of present data sample
%
%   Outputs:
%   dirCurrentName - directory only
%
%   Usage: 
%   name = generateName('csv', idx)

if (nargin < 3)
    idx = [];
end

name = [];

switch namingCase
    %     case 'matfilein'
    %         name = fullfile(getSetting('systemdir'), 'in.mat');
    %
    %     case 'matfilein-v7.3'
    %         name = fullfile(getSetting('systemdir'), 'in-v73.mat');
    %
    %     case 'matfileout'
    %         name = fullfile(getSetting('savedir'), getSetting('action'), 'out.mat');

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
        csv = generateName('csv', idx);
        time = strrep(idx.T, ':', ' ');
        time = strrep(time, '.', ',');
        un = strcat('(', num2str(idx.MsiID), ')');
        sample = idx.Sample;
        name = strjoin({csv, time, sample, un}, '_');

    case 'read'
        currentName = generateName('current', idx);
        name = fullfile(getSetting('savedir'), 'Cropped', currentName);

    case 'plot'
        currentName = generateName('current', idx);
        name = fullfile(getSetting('savedir'), getSetting('action'), currentName);

    case 'input'
        if isempty(getSetting('input'))
            name = fullfile('..', '..', '..', 'input');
        else
            name = getSetting('systemdir');
        end

    case 'output'
        if isempty(getSetting('output'))
            name = fullfile('..', '..', '..', 'output');
        else
            name = fullfile(getSetting('savedir'), getSetting('action'));
        end

    case 'id'
        name = fullfile(getSetting('systemdir'), 'ID.mat');

    case 'system'
        name = fullfile(getSetting('systemdir'), 'system.mat');

    otherwise
        name = fullfile(getSetting('savedir'), getSetting('action'), namingCase);

end

directory = fileparts(name);
mkNewDir(directory);

end

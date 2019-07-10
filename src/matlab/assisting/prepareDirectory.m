function [data, ID] = prepareDirectory(directory)
%PREPAREDIRECTORY Checks the data before reading, that everything is ok
%   if no argument, it runs in the current directory
%   else it runs in the designated directory
%
%     {'CoeffRef'   }
%     {'IMG'        }
%     {'IsCut'      }
%     {'IsFixed'    }
%     {'IsBenign'   }
%     {'Originx'    }
%     {'Originy'    }
%     {'Sample'     }
%     {'T'          }
%     {'Type'       }
%     {'MsiID'}
%     {'SpectrumFile'      }
%     {'data'       }
%     {'rmse'       }

%% setup
if nargin < 1
    directory = pwd;
end

oldDirectory = pwd;
if ~strcmp(directory, oldDirectory)
    cd(directory);
end

%% create data.mat
filetype = 'tif';
D = rdir(['*\*\*\*.', filetype]);
dataInfo = {D.name}';
data = struct('File', {}, 'Sample', {}, 'Type', {}, 'Camera', {}, 'Filename', {}, 'Extension', {}, 'MsiID', [], 'Count', [], 'HasPointImg', [], 'HasSpectralData', []);
count = 0;
for i = 1:length(dataInfo)
    dataSplit = strsplit(dataInfo{i}, {'\', '.'});
    data(i).File = dataInfo{i};
    data(i).Sample = dataSplit(1);
    data(i).Type = dataSplit(2);
    data(i).Camera = dataSplit(3);
    data(i).Filename = dataSplit(4);
    data(i).Extension = dataSplit(5);
    data(i).HasPointImg = false;
    data(i).HasSpectralData = false;
    if (i > 1) && ~(strcmp(data(i).Sample, data(i-1).Sample) && strcmp(data(i).Type, data(i-1).Type) && strcmp(data(i).Camera, data(i-1).Camera))
        count = 0;
    end
    count = count + 1;
    data(i).Count = ceil(count/8);
    data(i).MsiID = ceil(i/8);
end

pointImg = rdir('*\*\IMG*.jpg');
for i = 1:length(pointImg)
    dataSplit = strsplit(pointImg(i).name, {'\', '.'});
    idx2 = find(strcmp([data.Sample], dataSplit(1)));
    idx1 = find(strcmp([data.Type], dataSplit(2)));
    idx = intersect(idx1, idx2);
    for k = idx
        data(k).HasPointImg = true;
    end
end

fprintf('Total number of images = %d, of which:\n', data(end).MsiID);
fprintf('      unfixed = %d\n', sum(strcmp([data.Type], 'unfixed'))/8);
fprintf('      fixed = %d\n', sum(strcmp([data.Type], 'fixed'))/8);
fprintf('      cut = %d\n', sum(strcmp([data.Type], 'cut'))/8);

dirs = strcat([data.Sample], '\', [data.Type], '\', [data.Camera]);
for i = 1:length(dirs)
    if mod(numel(dir([dirs{i}, '\*.', filetype])), 8) > 0
        warning('MS subimages missing.')
        dirs(i)
        return
    end
end

%% create ID.mat
ID = struct('MsiID', [], 'Originx', [], 'Originy', [], 'IsBenign', [], 'T', {}, 'IsFixed', [], 'SpectrumFile', {}, 'IMG', {}, 'CoeffID', {});

csvList = rdir('*\*\*.csv');
if ~isempty(csvList)
    l = 1;
    for i = 1:length(csvList)
        dataSplit = strsplit(csvList(i).name, {'\', '.'});
        idx2 = find(strcmp([data.Sample], dataSplit(1)));
        idx1 = find(strcmp([data.Type], dataSplit(2)));
        idx = intersect(idx1, idx2);
        for k = idx
            data(k).HasSpectralData = true;
        end
        for k = idx(1:8:length(idx))
            [~, times, ~] = xlsread(csvList(i).name, 1, 'A3');
            times = strsplit(char(times), ',');
            if contains(csvList(i).name, {'notnormal', 'not_normal', 'cancer', 'shuyou'}, 'IgnoreCase', true)
                IsBenign = false;
            else
                IsBenign = true;
            end
            for n = 2:numel(times)
                ID(l).MsiID = data(k).MsiID;
                ID(l).Originx = 1;
                ID(l).Originy = 1;
                ID(l).SpectrumFile = csvList(i).name;
                if contains(csvList(i).name, {'notnormal', 'not_normal', 'cancer', 'shuyou'}, 'IgnoreCase', true)
                    ID(l).IsBenign = false;
                else
                    ID(l).IsBenign = true;
                end
                ID(l).IsBenign = IsBenign;
                ID(l).T = times{n};
                ID(l).Sample = data(ID(l).RgbID).Sample;
                ID(l).Type = data(ID(l).RgbID).Type;
                ID(l).IsFixed = ~contains(data(ID(l).RgbID).File, 'unfixed');
                ID(l).IsCut = contains(data(ID(l).RgbID).File, 'cut');
                l = l + 1;
            end
        end
    end
    a = {ID.MsiID};
    b = {ID.T};
    d = {ID.SpectrumFile};
    c = cellfun(@(x, y, z) [x, y, z], a', b', d', 'un', 0);
    [~, ii] = unique(c, 'stable');
    ID = ID(ii);
end


samples = unique([ID.Sample]);
for i = 1 : length(samples)
    for j = 1:length(ID)
        if strcmp(samples(i), ID(j).Sample)
            ID(j).SampleID = i;
        end
    end
end

measuredmat = dir('*MeasuredReflectance.mat');
if ~isempty(measuredmat)
    load(measuredmat.name);
    for i = 1:8:length(data)
        % sample name exists in the measument list
        if ~isempty(strfind({measuredReflectance.sampleName}, data(i).Sample))
            l = ceil(i/8);
            ID(l).MsiID = data(i).MsiID;
            ID(l).Originx = 1200;
            ID(l).Originy = 800;
            ID(l).IsFixed = false;
            ID(l).IsBenign = true;
            ID(l).Data = find([data.MsiID] == ID(k).MsiID);
        end
    end
end

%% finishing up
disp('Raw data check finished!')
cd(oldDirectory);

splits = strsplit(directory, '\');
if ~isempty(data)
    save(char(strcat('..\MATLAB\Data\', splits(end), '\data.mat')), 'data');
end
if ~isempty(ID)
    save(char(strcat('..\MATLAB\Data\', splits(end), '\ID.mat')), 'ID');
end

ID = orderfields(ID);
%writetable(struct2table(data), 'data.xlsx')
%ID  = table2struct( readtable('ID.xlsx') )
%     for i = 1:length(ID)
% ID(i).data = [ ID(i).data_1, ID(i).data_2,ID(i).data_3,ID(i).data_4,ID(i).data_5,ID(i).data_6,ID(i).data_7,ID(i).data_8];
% end

end

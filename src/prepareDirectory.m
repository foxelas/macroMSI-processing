function [ data, ID ] = prepareDirectory( directory )
%PREPAREDIRECTORY Checks the data before reading, that everything is ok
%   if no argument, it runs in the current directory
%   else it runs in the designated directory
% 
%     {'CoeffRef'   }
%     {'IMG'        }
%     {'IsCut'      }
%     {'IsFixed'    }
%     {'IsNormal'   }
%     {'Originx'    }
%     {'Originy'    }
%     {'Sample'     }
%     {'T'          }
%     {'Type'       }
%     {'UniqueCount'}
%     {'csvid'      }
%     {'data'       }
%     {'rmse'       }

%% setup
    if nargin < 1
        directory  = pwd;
    end
    
    oldDirectory = pwd;
    if ~strcmp(directory, oldDirectory)
        cd(directory);
    end    
    
%% create data.mat 
    filetype = 'tif';
    D = rdir(['*\*\*\*.', filetype]);
    dataInfo = {D.name}';
    data = struct('File', {}, 'Sample', {}, 'Type', {}, 'Camera', {}, 'Filename', {}, 'Extension', {}, 'UniqueCount', [], 'Count', [], 'HasPointImg', [], 'HasSpectralData', []);
    count = 0;
    for i=1:length(dataInfo)
        dataSplit = strsplit(dataInfo{i}, {'\', '.'});
        data(i).File = dataInfo{i};
        data(i).Sample = dataSplit(1);
        data(i).Type = dataSplit(2);
        data(i).Camera = dataSplit(3);
        data(i).Filename = dataSplit(4);
        data(i).Extension = dataSplit(5);
        data(i).HasPointImg = false;
        data(i).HasSpectralData = false;
        if (i > 1) && ~(strcmp(data(i).Sample,data(i-1).Sample) && strcmp(data(i).Type, data(i-1).Type) && strcmp(data(i).Camera,data(i-1).Camera))
            count = 0;
        end
        count = count + 1;
        data(i).Count = ceil(count/8);
        data(i).UniqueCount = ceil(i/8);
    end

    pointImg = rdir('*\*\IMG*.jpg');
    for i=1:length(pointImg)
        dataSplit = strsplit(pointImg(i).name, {'\', '.'});
        idx2 = find(strcmp([data.Sample], dataSplit(1)));
        idx1 = find(strcmp([data.Type], dataSplit(2)));
        idx = intersect(idx1,idx2);
        for k = idx
            data(k).HasPointImg = true;
        end
    end
    
    fprintf('Total number of images = %d, of which:\n', data(end).UniqueCount);
    fprintf('      unfixed = %d\n', sum(strcmp([data.Type], 'unfixed'))/8);
    fprintf('      fixed = %d\n', sum(strcmp([data.Type], 'fixed'))/8);
    fprintf('      cut = %d\n', sum(strcmp([data.Type], 'cut'))/8);

    dirs = strcat([data.Sample], '\', [data.Type], '\', [data.Camera]);
    for i=1:length(dirs)
        if mod(numel(dir( [dirs{i}, '\*.', filetype])),8) > 0 
            warning('MS subimages missing.')
            dirs(i)
            return
        end
    end
    
%% create ID.mat
    ID = struct('UniqueCount', [], 'Originx', [], 'Originy', [], 'IsNormal', [], 'T', {} , 'IsFixed', [], 'csvid', {}, 'IMG', {}, 'CoeffRef', {}, 'data', []);
        
    csvList = rdir('*\*\*.csv');
    if  ~isempty(csvList)
        l = 1;
        for i=1:length(csvList)
            dataSplit = strsplit(csvList(i).name, {'\', '.'});
            idx2 = find(strcmp([data.Sample], dataSplit(1)));
            idx1 = find(strcmp([data.Type], dataSplit(2)));
            idx = intersect(idx1,idx2);
            for k = idx
                data(k).HasSpectralData = true;
            end
            for k = idx(1:8:length(idx))
                [~, times,~] = xlsread(csvList(i).name, 1, 'A3');      
                times = strsplit(char(times), ',');
                if contains(csvList(i).name,{'notnormal', 'not_normal', 'cancer', 'shuyou'},'IgnoreCase',true)
                    isNormal = false;
                else 
                    isNormal = true;
                end
                for n = 2:numel(times)
                    ID(l).UniqueCount = data(k).UniqueCount;
                    ID(l).Originx = 1;
                    ID(l).Originy = 1;
                    ID(l).Csvid = csvList(i).name;
                    if contains(csvList(i).name,{'notnormal', 'not_normal', 'cancer', 'shuyou'},'IgnoreCase',true)
                        ID(l).IsNormal = false;
                    else 
                        ID(l).IsNormal = true;
                    end
                    ID(l).IsNormal = isNormal;
                    ID(l).T = times{n} ;
                    ID(l).Data = find([data.UniqueCount] == ID(l).UniqueCount);
                    ID(l).Sample = data(ID(l).Representative).Sample;
                    ID(l).Type = data(ID(l).Representative).Type;
                    ID(l).IsFixed = ~contains( data(ID(l).Representative).File, 'unfixed');
                    ID(l).IsCut = contains( data(ID(l).Representative).File, 'cut');
                    l = l + 1;
                end
            end
        end
        a= {ID.UniqueCount};
        b= {ID.T};
        d = {ID.Csvid};
        c=cellfun(@(x,y,z) [x y z],a', b', d', 'un',0);
        [~,ii]=unique(c,'stable');
        ID = ID(ii);
    end 
    
    measuredmat = dir('*MeasuredReflectance.mat');
    if ~isempty(measuredmat)
        load(measuredmat.name);
        for i=1:8:length(data)
            % sample name exists in the measument list 
            if ~isempty(strfind({measuredReflectance.sampleName}, data(i).Sample))
                l = ceil(i/8);
                ID(l).UniqueCount = data(i).UniqueCount;
                ID(l).Originx = 1200;
                ID(l).Originy = 800; 
                ID(l).IsFixed = false;
                ID(l).IsNormal = true;
                ID(l).Data = find([data.UniqueCount] == ID(k).UniqueCount);
            end
        end
    end
    %% finishing up
    disp('Raw data check finished!')
    cd(oldDirectory);   
    
    splits = strsplit(directory, '\');
    if ~isempty(data)
        save( char(strcat( '..\MATLAB\Data\', splits(end), '\data.mat')), 'data');
    end
    if ~isempty(ID)
        save( char(strcat( '..\MATLAB\Data\', splits(end), '\ID.mat')), 'ID'); 
    end
    
    ID = orderfields(ID);
    %writetable(struct2table(data), 'data.xlsx')
    %ID  = table2struct( readtable('ID.xlsx') )
%     for i = 1:length(ID)
% ID(i).data = [ ID(i).data_1, ID(i).data_2,ID(i).data_3,ID(i).data_4,ID(i).data_5,ID(i).data_6,ID(i).data_7,ID(i).data_8];
% end

end


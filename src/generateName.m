function [ savefilename, sampleName ] = generateName(options, filenameType, datax, idx, overrideName)
%GENERATEOUTNAME generates the name of the file,
%where the result will be saved.
%savefilename: the saveas filename icluding directory in the form (..\..\..\savename.extension)
%sampleName: the sample name 

    if (nargin < 2 )
        filenameType = 'default';
    end
    
    if (nargin < 3)
        datax = [];
        idx = [];
    end
    if (~isempty(datax) && isnumeric(datax.Sample))
        datax.Sample = num2str( datax.Sample);
    end
    sampleName = [];
    
    
    switch filenameType
        case 'default'
             extra = {};
             if strcmp(options.action, 'ReflectanceEstimation')
                 if strcmp(options.smoothingMatrixMethod, 'markovian' )
                    extra{end+1} = num2str(options.rho);
                 end 
                 if (options.withNoise) 
                     extra{end+1} = 'Noise';
                 else
                     extra{end+1} = 'NoNoise';
                 end
             end
             savefilename = strjoin([options.action, options.smoothingMatrixMethod, extra], '_');
             savefilename = fullfile(options.saveOptions.savedir, savefilename); 
             
        case 'matfilein'
            savefilename = fullfile(options.systemdir, 'in.mat');
            
        case 'matfileout'
             savefilename = fullfile(options.saveOptions.savedir, options.action, 'out.mat'); 
             
        case 'imagefile'
            if isempty(idx.T)
                savefilename = strjoin( [ datax.Camera{1}, num2str(idx.UniqueCount)] ,'_');
            else
                savefilename = strjoin( { datax.Camera{1}, idx.T, num2str(idx.UniqueCount)} ,'_');
            end
            savefilename = strrep(savefilename, ':', '-');
            savefilename = strrep(savefilename, '.', '-');
            
        case 'csvfile'
            if isempty(idx.Csvid)
                savefilename = num2str(idx.UniqueCount);
            else
                savefilename = strrep(idx.Csvid, '.csv', '');
                savefilename = strrep(savefilename, '\', '_'); 
            end
            
        case 'sample'
            splits = strsplit(idx.Csvid, '\');            
            savefilename = splits{1};
                  
        case 'plot+save'
            sampleName = generateName(options, 'imagefile',  datax, idx); 
            if (options.saveOptions.saveImages)
                fn = fullfile(options.saveOptions.savedir, options.action);
                savefilename = fullfile(fn, [options.action,'_', sampleName]);
            end
            
        case 'override'
            savefilename = fullfile(options.saveOptions.savedir, options.action, overrideName); 

        otherwise 
            warning('Could not generate name.')      
    end
       
end


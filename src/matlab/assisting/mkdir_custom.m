function [filepath] = mkdir_custom(filepath)
    filedir = fileparts(filepath);
    if ~exist(filedir, 'dir')
        mkdir(filedir);
        addpath(filedir);
    end
end
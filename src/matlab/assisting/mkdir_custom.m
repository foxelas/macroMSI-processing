function [] = mkdir_custom(filepath)
    if ~exist(filepath, 'dir')
        mkdir(filepath);
        addpath(filepath);
    end
end
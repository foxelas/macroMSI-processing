function [filepath] = mkNewDir(filepath)
%     MKNEWDIR creates a new directory 
% 
%     Usage:
%     [filepath] = mkNewDir(filepath)
 
filedir = fileparts(filepath);
if ~exist(filedir, 'dir')
    mkdir(filedir);
    addpath(filedir);
end
end
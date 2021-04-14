function [filepath] = mkNewDir(varargin)
%     MKNEWDIR creates a new directory
%
%     Usage:
%     [filepath] = mkNewDir(filepath)

if nargin == 1
    filepath = varargin{1};
else
    filepath = fullfile(varargin{:});
end
filedir = fileparts(filepath);
if ~exist(filedir, 'dir')
    mkdir(filedir);
    if  ~contains(filedir, getSetting('savedir'))
        addpath(filedir);
    end
end
end
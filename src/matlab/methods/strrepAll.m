function [outname] = strrepAll(inname)
%     STRREPALL fomats an inname to outname 
% 
%     Usage: 
%     [outname] = strrepAll(inname)

[~, outname] = fileparts(inname);
outname = strrep(outname, '\', ' ');
outname = strrep(outname, '_', ' ');
outname = strrep(outname, '.csv', '');
outname = strrep(outname, '.mat', '');

end
function [measuredSpectrum] = readSpectrum( csvpath, time )
%readSpectrum Reads the contents of a csvfile
%  and returns the measured spectrum at time 'time'

    if (nargin < 2)
        time = 'no time';
    end
    
    outstruct = delimread(csvpath,',','raw');
    
    % if time is specified, find the respective column, else it's the 5th
    idx = max( [5, find(not(cellfun('isempty',  strfind( outstruct.raw(3,:),  time) )))] );

    measuredSpectrum =  cell2mat(cellfun(@str2num,  outstruct.raw(26:end, idx),'un',0));

end


function spectrum = readSpectrum(csvpath, time)
%     READSPECTRUM Read the contents of a csvfile to a vector of spectral data
%
%     Inputs:
%     csvpath - the path and filename of the csvfile containing multiple
%     spectral info
%     time - the time value for the spectal info
%
%     Outputs:
%     spectrum - the read spectrum
%
%     Usage:
%     spectrum = readSpectrum(csvpath, time)
%     spectrum = readSpectrum(csvpath)


if (nargin < 2)
    time = 'no time';
end

outstruct = delimread(csvpath, ',', 'raw');

% if time is specified, find the respective column, else it's the 5th
idx = max([5, find(not(cellfun(@(x) isempty(x), strfind(outstruct.raw(3, :), time))))]);

spectrum = cell2mat(cellfun(@str2num, outstruct.raw(26:end, idx), 'un', 0));

end

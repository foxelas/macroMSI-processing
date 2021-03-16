function [values, valueNames, additionalValues] = getExpectedValues(name, option)
%GETEXPECTEDVALUES fetches expected values
% 
%   [values, valueNames, additionalValues] = getExpectedValues() returns 
%   the expected values for colorchart spectra
%
%   [values] = getExpectedValues('colorchartRGB') returns the expected
%   values for colorchart RGB space
%
%   [values] = getExpectedValues('colorchartLab') returns the
%   expected values for colorchart Lab space
%
%   [values] = getExpectedValues('colorchartOrder', 
%   'capture_average_comparison') returns the expected values for
%   colorchart patch order 
%

    if nargin < 1
        name = 'colorchartSpectra';
    end 

    valueNames = [];
    additionalValues = [];
    switch name 
        case 'colorchartSpectra'
            filename = fullfile(getSetting('systemdir'), 'ColorChecker_RGB_and_spectra.txt');
            outstruct = delimread(filename, '\t', {'text', 'num'});
            valueNames = outstruct.text;
            valueNames = valueNames(2:length(valueNames));
            additionalValues = outstruct.num(1, :);
            values = outstruct.num(2:end, :);

        case 'colorchartRGB'
            filename = fullfile(getSetting('systemdir'), 'ColorCheckerMicro_Matte_RGB_values.txt');
            outstruct = delimread(filename, '\t', 'num');
            values = outstruct.num;
        
        case 'colorchartLab'
            filename = fullfile(getSetting('systemdir'), 'ColorCheckerMicro_Matte_Lab_values.txt');
            outstruct = delimread(filename, '\t', 'num');
            values = outstruct.num;
            
        case 'colorchartOrder'
            configuration = option; 
            if ~strcmp(configuration, 'singleLightFar') && ~strcmp(configuration, 'testStomach')
                configuration = 'colorchart';
            end 
            outstruct = delimread(fullfile(getSetting('datasetSettingsDir'), strcat(configuration, 'PatchOrder.txt')), '\t', 'text');
            values = outstruct.text;

        otherwise 
            error('Unsupported name.')
    end 
    
end 
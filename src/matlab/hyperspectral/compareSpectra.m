function [T, adjusted, alphaCoeff] = compareSpectra(expected, measured, lineNames)
%COMPARESPECTRA compares expected and measured spectra
%
%   T = compareSpectra(expected, measured, lineNames) returns a table with
%   comparison differences for three metrics, Goodness-Of-Fit, Normalized
%   Mean Square Error and Root Mean Square Error
%
%   [T, adjusted, alpha] = compareSpectra(expected, measured, lineNames) 
%   returns table T and also adjusted values after adjustment
%

%% Standard (expected) color patch spectra
plots(4, @plotColorChartSpectra, expected, lineNames, 'expected');
plots(5, @plotColorChartSpectra, measured, lineNames, 'measured');

m = size(measured, 2);
x = getWavelengths(m, 'babel');
expected = expected(:, x);

% Limit to Range [420,730]nm 
measured = measured(:,5:end);
expected = expected(:,5:end);

difference = expected - measured;
plots(6, @plotColorChartSpectra, difference, lineNames, 'difference');

gofs = applyRowFunc(@goodnessOfFit, measured, expected);
nmses = applyRowFunc(@nmse, measured, expected);
rmses = applyRowFunc(@rmse, measured, expected);

%% Standard (expected) color patch spectra after adjustment
[adjusted, alphaCoeff] = adjustSpectra(measured, lineNames, 'toRatio');

plots(7, @plotColorChartSpectra, adjusted, lineNames, 'measured-adjusted');

differenceAdjusted = expected - adjusted;
plots(8, @plotColorChartSpectra, differenceAdjusted, lineNames, 'difference-adjusted');

gofsAdj = applyRowFunc(@goodnessOfFit, adjusted, expected);
nmsesAdj = applyRowFunc(@nmse, adjusted, expected);
rmsesAdj = applyRowFunc(@rmse, adjusted, expected);

T = table(lineNames, gofs, nmses, rmses, gofsAdj, nmsesAdj, rmsesAdj, ...
    'VariableNames', {'Patch', 'GoF', 'NMSE', 'RMSE', 'AdjGoF', 'AdjNMSE', 'AdjRMSE'});

end
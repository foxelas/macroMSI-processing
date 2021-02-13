function [T, adjusted] = compareSpectra(expected, measured, lineNames)
%COMPARESPECTRA compares expected and measured spectra 
%
%   T = compareSpectra(expected, measured, lineNames) returns a table with
%   comparison differences for three metrics, Goodness-Of-Fit, Normalized
%   Mean Square Error and Root Mean Square Error 
%
%   [T, adjusted] = compareSpectra(expected, measured, lineNames) returns
%   table T and also adjusted values after adjustment 
%

%% Standard (expected) color patch spectra
plots(4, @plotColorChartSpectra, expected, lineNames, 'expected');
plots(5, @plotColorChartSpectra, measured, lineNames, 'measured');

m = size(measured, 2);
x = getWavelengths(m, 'index');
expected = expected(:,x);

difference = expected - measured;
plots(6, @plotColorChartSpectra, difference, lineNames, 'difference');

gofs = applyRowFunc(measured, expected, @goodnessOfFit);
nmses = applyRowFunc(measured, expected, @nmse);
rmses = applyRowFunc(measured, expected, @rmse);

%% Standard (expected) color patch spectra after adjustment
[adjusted, ~] = adjustSpectra(measured, lineNames);

plots(7, @plotColorChartSpectra, adjusted, lineNames, 'measured-adjusted');

differenceAdjusted = expected - adjusted;
plots(8, @plotColorChartSpectra, differenceAdjusted, lineNames, 'difference-adjusted');

gofsAdj = applyRowFunc(adjusted, expected, @goodnessOfFit);
nmsesAdj = applyRowFunc(adjusted, expected, @nmse);
rmsesAdj = applyRowFunc(adjusted, expected, @rmse);

T = table(lineNames, gofs',nmses',rmses',gofsAdj',nmsesAdj',rmsesAdj', ...
    'VariableNames' , {'Patch', 'GoF', 'NMSE', 'RMSE', 'AdjGoF', 'AdjNMSE', 'AdjRMSE'});

end
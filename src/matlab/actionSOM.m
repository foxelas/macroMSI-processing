outfile = matfile(fullfile(options.saveOptions.savedir, outputFolderMap('features'), 'out.mat'));
color = outfile.EstimatedSpectra';
mmlbp = outfile.MMLbpFeatures; 
texture = [ mmlbp{1},  mmlbp{2} ]';

%% only color 
inputs  = color;
dimension1 = 10;
dimension2 = 10;
getSom(inputs, dimension1, dimension2, options.saveOptions, 'spect');

%% color and texture 
inputs = [ color; texture ]; 

% Create a Self-Organizing Map
dimension1 = 10;
dimension2 = 10;
getSom(inputs, dimension1, dimension2, options.saveOptions, 'spect+mmlbp');
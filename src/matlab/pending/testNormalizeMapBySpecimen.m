% files = {'F:\temp\mspi\matfiles\macroMSI\24-Nov-2020_Consistent_with_Veruca_Vulgaris.mat', ...
%     'F:\temp\mspi\matfiles\macroMSI\23-Nov-2020_Spitz_Nevus.mat', ...
%     'F:\temp\mspi\matfiles\macroMSI\23-Nov-2020_Melanocytic_Nevus.mat', ...
%     'F:\temp\mspi\matfiles\macroMSI\24-Nov-2020_Basal_Cell_Carcinoma.mat', ...
%     'F:\temp\mspi\matfiles\macroMSI\24-Nov-2020_Pigmented_Seborrheic_Keratosis.mat'};
 

close all; 
clear all; 
load('F:\temp\mspi\matfiles\macroMSI\24-Nov-2020_Consistent_with_Veruca_Vulgaris.mat');
setSetting('savedir',  mkNewDir('..\..\..\output\', 'spieNormalizeBySpecimen\'));
actionGetMaps;

close all; 
clear all; 
load('F:\temp\mspi\matfiles\macroMSI\23-Nov-2020_Spitz_Nevus.mat');
setSetting('savedir',  mkNewDir('..\..\..\output\', 'spieNormalizeBySpecimen\'));
actionGetMaps;

close all; 
clear all; 
load('F:\temp\mspi\matfiles\macroMSI\23-Nov-2020_Melanocytic_Nevus.mat');
setSetting('savedir',  mkNewDir('..\..\..\output\', 'spieNormalizeBySpecimen\'));
actionGetMaps;

close all; 
clear all; 
load('F:\temp\mspi\matfiles\macroMSI\24-Nov-2020_Basal_Cell_Carcinoma.mat');
setSetting('savedir',  mkNewDir('..\..\..\output\', 'spieNormalizeBySpecimen\'));
actionGetMaps;

close all; 
clear all; 
load('F:\temp\mspi\matfiles\macroMSI\24-Nov-2020_Pigmented_Seborrheic_Keratosis.mat');
setSetting('savedir',  mkNewDir('..\..\..\output\', 'spieNormalizeBySpecimen\'));
actionGetMaps;

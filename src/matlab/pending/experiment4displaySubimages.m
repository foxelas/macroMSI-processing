close all; clc; 

k = 20;
msiType = 'extended'; %'extended'; % 'max';
 
% for all pixels without bg removal 
removebg = true; 
[msi, ~, ~, ~, ~, ~] = getImage(k, options, msiType, removebg, false);
plotMSI(msi, 1, options.saveOptions);
title('Msi (Extended)')
[difMsi, avgMsi] = getDifMSI(msi, 'toAverage') ;
plotMSI(avgMsi, 2, options.saveOptions);
title('Average')
plotMSI(difMsi, 3, options.saveOptions);
title('Difference to Average')
[difMsi, ~] = getDifMSI(msi, 'toNextBand') ;
plotMSI(difMsi, 4, options.saveOptions);
title('Difference to Next Band')
[normMsi] = getNormMSI(msi);
plotMSI(normMsi, 5, options.saveOptions);
title('Normalized Msi') 
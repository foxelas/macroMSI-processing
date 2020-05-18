close all; clc; 

k = 20;
msiType = 'extended'; %'extended'; % 'max';
 
% for all pixels without bg removal 
removebg = true; 
[msi, ~, ~, ~, ~, ~] = getImage(k, msiType, removebg, false);
plotMSI(msi, 1);
title('Msi (Extended)')
[difMsi, avgMsi] = getDifMSI(msi, 'toAverage') ;
plotMSI(avgMsi, 2);
title('Average')
plotMSI(difMsi, 3);
title('Difference to Average')
[difMsi, ~] = getDifMSI(msi, 'toNextBand') ;
plotMSI(difMsi, 4);
title('Difference to Next Band')
[normMsi] = getNormMSI(msi);
plotMSI(normMsi, 5);
title('Normalized Msi') 
close all; clc;

k = 20;
msiType = 'extended'; %'extended'; % 'max';

% for all pixels without bg removal
removebg = true;
[msi, ~, ~, ~, ~, ~] = getImage(k, msiType, removebg, false);
plotFunWrapper(1, @plotMSI, msi);
title('Msi (Extended)')
[difMsi, avgMsi] = getDifMSI(msi, 'toAverage');
plotFunWrapper(2, @plotMSI, avgMsi);
title('Average')
plotFunWrapper(3, @plotMSI, difMsi);
title('Difference to Average')
[difMsi, ~] = getDifMSI(msi, 'toNextBand');
plotFunWrapper(4, @plotMSI, difMsi);
title('Difference to Next Band')
[normMsi] = getNormMSI(msi);
plotFunWrapper(5, @plotMSI, normMsi);
title('Normalized Msi')
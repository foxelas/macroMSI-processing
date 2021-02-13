close all; clc;

k = 20;
msiType = 'extended'; %'extended'; % 'max';

% for all pixels without bg removal
removebg = true;
[msi, ~, ~, ~, ~, ~] = getImage(k, msiType, removebg, false);
plots(1, @plotMSI, msi);
title('Msi (Extended)')
[difMsi, avgMsi] = getDifMSI(msi, 'toAverage');
plots(2, @plotMSI, avgMsi);
title('Average')
plots(3, @plotMSI, difMsi);
title('Difference to Average')
[difMsi, ~] = getDifMSI(msi, 'toNextBand');
plots(4, @plotMSI, difMsi);
title('Difference to Next Band')
[normMsi] = getNormMSI(msi);
plots(5, @plotMSI, normMsi);
title('Normalized Msi')
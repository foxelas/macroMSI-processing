xlsxfilename = 'hsiEvaluation.xlsx';

load('D:\temp\Google Drive\titech\research\experiments\output\hsi\singleLightFar_byAverage\last_run.mat')
white95Idx = find(strcmp(spectraColorOrder,  'white 9.5 (.05 D)')); 
white95Val = 0.8;
a = mean(reorderedSpectralVals(white95Idx,(end-20):end)) / white95Val;
fprintf('Values adjusted so that white 9.5 (.05 D) line is assinged to value 0.8 \nwith division by alpha = %.3f \n', a);
adjustedReorderedSpectralVals = reorderedSpectralVals / a;
gofs2 = applyFuncOnRows(adjustedReorderedSpectralVals, expectedSpectra, @goodnessOfFit);
nmses2 = applyFuncOnRows(adjustedReorderedSpectralVals, expectedSpectra, @nmse);
lessVals = adjustedReorderedSpectralVals(:, 4:end );
expectedSpectra = expectedSpectra(:,4:end);
gofs3 = applyFuncOnRows(lessVals, expectedSpectra, @goodnessOfFit);
nmses3 = applyFuncOnRows(lessVals, expectedSpectra, @nmse);

writematrix('Conf1Norm1',xlsxfilename,'Sheet',1,'Range','A2');
writematrix('Conf1Norm1',xlsxfilename,'Sheet',1,'Range','A12');
writematrix(nmses,xlsxfilename,'Sheet',1,'Range','B2:Y2');
writematrix(gofs,xlsxfilename,'Sheet',1,'Range','B12:Y12');
writematrix('Conf1Norm1',xlsxfilename,'Sheet',2,'Range','A2');
writematrix('Conf1Norm1',xlsxfilename,'Sheet',2,'Range','A12');
writematrix(nmses2,xlsxfilename,'Sheet',2,'Range','B2:Y2');
writematrix(gofs2,xlsxfilename,'Sheet',2,'Range','B12:Y12');
writematrix('Conf1Norm1',xlsxfilename,'Sheet',3,'Range','A2');
writematrix('Conf1Norm1',xlsxfilename,'Sheet',3,'Range','A12');
writematrix(nmses3,xlsxfilename,'Sheet',3,'Range','B2:Y2');
writematrix(gofs3,xlsxfilename,'Sheet',3,'Range','B12:Y12');

load('D:\temp\Google Drive\titech\research\experiments\output\hsi\singleLightFar_byPixel_noSmoothing\last_run.mat')
white95Idx = find(strcmp(spectraColorOrder,  'white 9.5 (.05 D)')); 
white95Val = 0.8;
a = mean(reorderedSpectralVals(white95Idx,(end-20):end)) / white95Val;
fprintf('Values adjusted so that white 9.5 (.05 D) line is assinged to value 0.8 \nwith division by alpha = %.3f \n', a);
adjustedReorderedSpectralVals = reorderedSpectralVals / a;
gofs2 = applyFuncOnRows(adjustedReorderedSpectralVals, expectedSpectra, @goodnessOfFit);
nmses2 = applyFuncOnRows(adjustedReorderedSpectralVals, expectedSpectra, @nmse);
lessVals = adjustedReorderedSpectralVals(:, 4:end );
expectedSpectra = expectedSpectra(:,4:end);
gofs3 = applyFuncOnRows(lessVals, expectedSpectra, @goodnessOfFit);
nmses3 = applyFuncOnRows(lessVals, expectedSpectra, @nmse);

writematrix('Conf1Norm2',xlsxfilename,'Sheet',1,'Range','A3');
writematrix('Conf1Norm2',xlsxfilename,'Sheet',1,'Range','A13');
writematrix(nmses,xlsxfilename,'Sheet',1,'Range','B3:Y3');
writematrix(gofs,xlsxfilename,'Sheet',1,'Range','B13:Y13');
writematrix('Conf1Norm2',xlsxfilename,'Sheet',2,'Range','A3');
writematrix('Conf1Norm2',xlsxfilename,'Sheet',2,'Range','A13');
writematrix(nmses2,xlsxfilename,'Sheet',2,'Range','B3:Y3');
writematrix(gofs2,xlsxfilename,'Sheet',2,'Range','B13:Y13');
writematrix('Conf1Norm2',xlsxfilename,'Sheet',3,'Range','A3');
writematrix('Conf1Norm2',xlsxfilename,'Sheet',3,'Range','A13');
writematrix(nmses3,xlsxfilename,'Sheet',3,'Range','B3:Y3');
writematrix(gofs3,xlsxfilename,'Sheet',3,'Range','B13:Y13');


load('D:\temp\Google Drive\titech\research\experiments\output\hsi\singleLightClose_byAverage\last_run.mat')
white95Idx = find(strcmp(spectraColorOrder,  'white 9.5 (.05 D)')); 
white95Val = 0.8;
a = mean(reorderedSpectralVals(white95Idx,(end-20):end)) / white95Val;
fprintf('Values adjusted so that white 9.5 (.05 D) line is assinged to value 0.8 \nwith division by alpha = %.3f \n', a);
adjustedReorderedSpectralVals = reorderedSpectralVals / a;
gofs2 = applyFuncOnRows(adjustedReorderedSpectralVals, expectedSpectra, @goodnessOfFit);
nmses2 = applyFuncOnRows(adjustedReorderedSpectralVals, expectedSpectra, @nmse);
lessVals = adjustedReorderedSpectralVals(:, 4:end );
expectedSpectra = expectedSpectra(:,4:end);
gofs3 = applyFuncOnRows(lessVals, expectedSpectra, @goodnessOfFit);
nmses3 = applyFuncOnRows(lessVals, expectedSpectra, @nmse);

writematrix('Conf2Norm1',xlsxfilename,'Sheet',1,'Range','A4');
writematrix('Conf2Norm1',xlsxfilename,'Sheet',1,'Range','A14');
writematrix(nmses,xlsxfilename,'Sheet',1,'Range','B4:Y4');
writematrix(gofs,xlsxfilename,'Sheet',1,'Range','B14:Y14');
writematrix('Conf2Norm1',xlsxfilename,'Sheet',2,'Range','A4');
writematrix('Conf2Norm1',xlsxfilename,'Sheet',2,'Range','A14');
writematrix(nmses2,xlsxfilename,'Sheet',2,'Range','B4:Y4');
writematrix(gofs2,xlsxfilename,'Sheet',2,'Range','B14:Y14');
writematrix('Conf2Norm1',xlsxfilename,'Sheet',3,'Range','A4');
writematrix('Conf2Norm1',xlsxfilename,'Sheet',3,'Range','A14');
writematrix(nmses3,xlsxfilename,'Sheet',3,'Range','B4:Y4');
writematrix(gofs3,xlsxfilename,'Sheet',3,'Range','B14:Y14');


load('D:\temp\Google Drive\titech\research\experiments\output\hsi\singleLightClose_byPixel_noSmoothing\last_run.mat')
white95Idx = find(strcmp(spectraColorOrder,  'white 9.5 (.05 D)')); 
white95Val = 0.8;
a = mean(reorderedSpectralVals(white95Idx,(end-20):end)) / white95Val;
fprintf('Values adjusted so that white 9.5 (.05 D) line is assinged to value 0.8 \nwith division by alpha = %.3f \n', a);
adjustedReorderedSpectralVals = reorderedSpectralVals / a;
gofs2 = applyFuncOnRows(adjustedReorderedSpectralVals, expectedSpectra, @goodnessOfFit);
nmses2 = applyFuncOnRows(adjustedReorderedSpectralVals, expectedSpectra, @nmse);
lessVals = adjustedReorderedSpectralVals(:, 4:end );
expectedSpectra = expectedSpectra(:,4:end);
gofs3 = applyFuncOnRows(lessVals, expectedSpectra, @goodnessOfFit);
nmses3 = applyFuncOnRows(lessVals, expectedSpectra, @nmse);

writematrix('Conf2Norm2',xlsxfilename,'Sheet',1,'Range','A5');
writematrix('Conf2Norm2',xlsxfilename,'Sheet',1,'Range','A15');
writematrix(nmses,xlsxfilename,'Sheet',1,'Range','B5:Y5');
writematrix(gofs,xlsxfilename,'Sheet',1,'Range','B15:Y15');
writematrix('Conf2Norm2',xlsxfilename,'Sheet',2,'Range','A5');
writematrix('Conf2Norm2',xlsxfilename,'Sheet',2,'Range','A15');
writematrix(nmses2,xlsxfilename,'Sheet',2,'Range','B5:Y5');
writematrix(gofs2,xlsxfilename,'Sheet',2,'Range','B15:Y15');
writematrix('Conf2Norm2',xlsxfilename,'Sheet',3,'Range','A5');
writematrix('Conf2Norm2',xlsxfilename,'Sheet',3,'Range','A15');
writematrix(nmses3,xlsxfilename,'Sheet',3,'Range','B5:Y5');
writematrix(gofs3,xlsxfilename,'Sheet',3,'Range','B15:Y15');

load('D:\temp\Google Drive\titech\research\experiments\output\hsi\doubleLightClose_byAverage\last_run.mat')
white95Idx = find(strcmp(spectraColorOrder,  'white 9.5 (.05 D)')); 
white95Val = 0.8;
a = mean(reorderedSpectralVals(white95Idx,(end-20):end)) / white95Val;
fprintf('Values adjusted so that white 9.5 (.05 D) line is assinged to value 0.8 \nwith division by alpha = %.3f \n', a);
adjustedReorderedSpectralVals = reorderedSpectralVals / a;
gofs2 = applyFuncOnRows(adjustedReorderedSpectralVals, expectedSpectra, @goodnessOfFit);
nmses2 = applyFuncOnRows(adjustedReorderedSpectralVals, expectedSpectra, @nmse);
lessVals = adjustedReorderedSpectralVals(:, 4:end );
expectedSpectra = expectedSpectra(:,4:end);
gofs3 = applyFuncOnRows(lessVals, expectedSpectra, @goodnessOfFit);
nmses3 = applyFuncOnRows(lessVals, expectedSpectra, @nmse);

writematrix('Conf3Norm1',xlsxfilename,'Sheet',1,'Range','A6');
writematrix('Conf3Norm1',xlsxfilename,'Sheet',1,'Range','A16');
writematrix(nmses,xlsxfilename,'Sheet',1,'Range','B6:Y6');
writematrix(gofs,xlsxfilename,'Sheet',1,'Range','B16:Y16');
writematrix('Conf3Norm1',xlsxfilename,'Sheet',2,'Range','A6');
writematrix('Conf3Norm1',xlsxfilename,'Sheet',2,'Range','A16');
writematrix(nmses2,xlsxfilename,'Sheet',2,'Range','B6:Y6');
writematrix(gofs2,xlsxfilename,'Sheet',2,'Range','B16:Y16');
writematrix('Conf3Norm1',xlsxfilename,'Sheet',3,'Range','A6');
writematrix('Conf3Norm1',xlsxfilename,'Sheet',3,'Range','A16');
writematrix(nmses3,xlsxfilename,'Sheet',3,'Range','B6:Y6');
writematrix(gofs3,xlsxfilename,'Sheet',3,'Range','B16:Y16');


load('D:\temp\Google Drive\titech\research\experiments\output\hsi\doubleLightClose_byPixel_noSmoothing\last_run.mat')
white95Idx = find(strcmp(spectraColorOrder,  'white 9.5 (.05 D)')); 
white95Val = 0.8;
a = mean(reorderedSpectralVals(white95Idx,(end-20):end)) / white95Val;
fprintf('Values adjusted so that white 9.5 (.05 D) line is assinged to value 0.8 \nwith division by alpha = %.3f \n', a);
adjustedReorderedSpectralVals = reorderedSpectralVals / a;
gofs2 = applyFuncOnRows(adjustedReorderedSpectralVals, expectedSpectra, @goodnessOfFit);
nmses2 = applyFuncOnRows(adjustedReorderedSpectralVals, expectedSpectra, @nmse);
lessVals = adjustedReorderedSpectralVals(:, 4:end );
expectedSpectra = expectedSpectra(:,4:end);
gofs3 = applyFuncOnRows(lessVals, expectedSpectra, @goodnessOfFit);
nmses3 = applyFuncOnRows(lessVals, expectedSpectra, @nmse);

writematrix('Conf3Norm2',xlsxfilename,'Sheet',1,'Range','A7');
writematrix('Conf3Norm2',xlsxfilename,'Sheet',1,'Range','A17');
writematrix(nmses,xlsxfilename,'Sheet',1,'Range','B7:Y7');
writematrix(gofs,xlsxfilename,'Sheet',1,'Range','B17:Y17');
writematrix('Conf3Norm2',xlsxfilename,'Sheet',2,'Range','A7');
writematrix('Conf3Norm2',xlsxfilename,'Sheet',2,'Range','A17');
writematrix(nmses2,xlsxfilename,'Sheet',2,'Range','B7:Y7');
writematrix(gofs2,xlsxfilename,'Sheet',2,'Range','B17:Y17');
writematrix('Conf3Norm2',xlsxfilename,'Sheet',3,'Range','A7');
writematrix('Conf3Norm2',xlsxfilename,'Sheet',3,'Range','A17');
writematrix(nmses3,xlsxfilename,'Sheet',3,'Range','B7:Y7');
writematrix(gofs3,xlsxfilename,'Sheet',3,'Range','B17:Y17');



